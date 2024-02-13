use crate::crypto::blake2b512::blake2b512;
use crate::field::field::Field;
use crate::field::field_element::FieldElement;
use crate::field::polynomial::Polynomial;
use crate::merkle_root::MerkleRoot;
use crate::proof_stream::PROOF_BYTES;
use crate::stark::{StarkProofStreamEnum, StarkProofStream};
use crate::utils::bit_iter::BitIter;
use crate::utils::bytes::Bytes;

pub struct FRI<'a> {
  omega: FieldElement<'a>,
  offset: FieldElement<'a>,
  field: &'a Field,
  domain_length: usize,
  expansion_factor: usize,
  num_colinearity_tests: usize,
}

impl<'a> FRI<'a> {
  pub fn new (
    offset: FieldElement<'a>,
    omega: FieldElement<'a>,
    domain_length: usize,
    expansion_factor: usize,
    num_colinearity_tests: usize,
  ) -> Self {
    Self {
      omega,
      offset,
      field: omega.field,
      domain_length,
      expansion_factor,
      num_colinearity_tests,
    }
  }

  pub fn num_rounds (&self) -> usize {
    let mut codeword_length = self.domain_length;
    let mut num_rounds = 0;

    while (codeword_length > self.expansion_factor) && (codeword_length > (4 * self.num_colinearity_tests)) {
      codeword_length /= 2;
      num_rounds += 1;
    }

    num_rounds
  }

  pub fn eval_domain (&self) -> Vec<FieldElement<'a>> {
    (0..self.domain_length)
      .map(|i| {
        self.offset * (self.omega.clone() ^ i)
      })
      .collect()
  }

  fn sample_index (&self, bytes: Bytes, size: usize) -> usize {
    assert_ne!(size, 0, "modulo zero is impossible");

    let iter: BitIter<usize> = size.into();
    let bit = iter.bit_index().unwrap(); // `None` only if zero, asserted above
    let bytes_num = bit / 8 + 1;

    let bytes = bytes.bytes();
    let len = bytes.len();

    let start_i = if bytes_num > len {
      0
    } else {
      len-bytes_num
    };
    let res = bytes[start_i..len]
      .iter()
      .fold(0, |acc, b| {
        let res = (acc << 8) ^ (*b as usize);
        res
      });

    res % size
  }

  pub fn sample_indices (
    &self,
    seed: Bytes,
    size: usize,
    reduced_size: usize,
    number: usize,
  ) -> Vec<usize> {
    assert!(number <= 2 * reduced_size, "Not enough entropy in indices with reference to last codeword");
    assert!(number <= reduced_size, "Cannot sample more indices than available in the last codeword");

    let mut indices = Vec::with_capacity(number);
    let mut reduced_indices = Vec::with_capacity(number);
    let mut counter = 0;

    while indices.len() < number {
      let bytes = seed.clone() + Bytes::new(vec![0; counter]);

      let index = self.sample_index(blake2b512(bytes), size);
      let reduced_index = index % reduced_size;
      counter += 1;

      if !reduced_indices.contains(&reduced_index) {
        indices.push(index);
        reduced_indices.push(reduced_index);
      }
    }

    indices
  }

  pub fn commit (
    &self,
    codeword: &Vec<FieldElement<'a>>,
    proof_stream: &mut StarkProofStream<'a>,
  ) -> Vec<Vec<FieldElement<'a>>> {
    let one = self.field.one();
    let two_inv = FieldElement::new(self.field, 2).inverse();
    let mut omega = self.omega;
    let mut offset = self.offset;

    let num_rounds = self.num_rounds();

    let mut codewords = Vec::with_capacity(num_rounds);
    let mut codeword = codeword.to_vec();

    for r in 0..num_rounds {
      /*
        TODO:
          serialization differs
          => merkle root is different
          => proof_stream is different
          => fiat_shamir_prover is different
          => last codeword is different
          BUT first codeword is the same. Maybe that's the problem?
      */

      // compute and send Merkle root
      let root = MerkleRoot::commit(&codeword);
      proof_stream.push(StarkProofStreamEnum::Root(root));

      // prepare next round, if necessary
      if r == num_rounds - 1 {
        break
      }

      // get challenge
      let challenge = proof_stream.fiat_shamir_prover(PROOF_BYTES);
      let alpha = self.field.sample(&challenge);

      codewords.push(codeword.clone());

      // split and fold
      let codeword_half_len = codeword.len() / 2;
      codeword = (0..codeword_half_len)
        .map(|i| {
          let alpha_by_offset = alpha / (offset * (omega ^ i));
          let first_half = (one + alpha_by_offset) * codeword[i].clone();
          let second_half = (one - alpha_by_offset) * codeword[codeword_half_len + i].clone();
          two_inv * (first_half + second_half)
        })
        .collect();

      omega = omega ^ 2_usize;
      offset = offset ^ 2_usize;
    }

    // send last codeword
    proof_stream.push(StarkProofStreamEnum::Codeword(codeword.clone()));

    // collect last codeword too
    codewords.push(codeword);

    codewords
  }

  pub fn query (
    &self,
    codeword_current: &[FieldElement<'a>],
    codeword_next: &[FieldElement<'a>],
    indices_c: &[usize],
    proof_stream: &mut StarkProofStream<'a>,
  ) -> Vec<usize> {
    // infer a and b indices
    let indices_a = indices_c.to_vec();
    let indices_b = indices_c
      .iter()
      .map(|i| i + codeword_current.len() / 2)
      .collect::<Vec<_>>();

    // reveal leafs
    for s in 0..self.num_colinearity_tests {
      proof_stream.push(StarkProofStreamEnum::Leafs((
        *codeword_current.get(*indices_a.get(s).unwrap()).unwrap(),
        *codeword_current.get(*indices_b.get(s).unwrap()).unwrap(),
        *codeword_next.get(*indices_c.get(s).unwrap()).unwrap(),
      )))
    }

    // reveal auth path
    for s in 0..self.num_colinearity_tests {
      proof_stream.push(StarkProofStreamEnum::Path(MerkleRoot::open(*indices_a.get(s).unwrap(), codeword_current)));
      proof_stream.push(StarkProofStreamEnum::Path(MerkleRoot::open(*indices_b.get(s).unwrap(), codeword_current)));
      proof_stream.push(StarkProofStreamEnum::Path(MerkleRoot::open(*indices_c.get(s).unwrap(), codeword_next)));
    }

    let mut indices_ab = indices_a;
    indices_ab.extend(indices_b);

    indices_ab
  }

  pub fn prove (
    &self,
    codeword: &Vec<FieldElement<'a>>,
    proof_stream: &mut StarkProofStream<'a>,
  ) -> Vec<usize> {
    assert_eq!(
      self.domain_length,
      codeword.len(),
      "Length of the domain doesnt match the length of initial codeword",
    );

    let codewords = self.commit(codeword, proof_stream);

    let top_level_indices = self.sample_indices(
      proof_stream.fiat_shamir_prover(PROOF_BYTES),
      codewords.get(1).unwrap().len(),
      codewords.last().unwrap().len(),
      self.num_colinearity_tests,
    );
    let mut indices = top_level_indices.clone();

    for i in 0..codewords.len()-1 {
      let el = codewords.get(i).unwrap();

      indices = indices
        .iter()
        .map(|i| i % (el.len() / 2))
        .collect();

      self.query(
        el,
        codewords.get(i + 1).unwrap(),
        &indices,
        proof_stream,
      );
    }

    top_level_indices
  }

  pub fn verify (
    &self,
    proof_stream: &mut StarkProofStream<'a>,
    polynomial_values: &mut Vec<(usize, FieldElement<'a>)>,
  ) -> Result<(), String> {
    let mut omega = self.omega;
    let mut offset = self.offset;

    let num_rounds = self.num_rounds();

    let mut roots = Vec::with_capacity(num_rounds);
    let mut alphas = Vec::with_capacity(num_rounds);

    for r in 0..num_rounds {
      let root = match proof_stream.pull() {
        Some(StarkProofStreamEnum::Root(root)) => root,
        _ => {
          return Err(format!("expected {r}th item in proof stream to be root"));
        }
      };

      roots.push(root);
      alphas.push(self.field.sample(&proof_stream.fiat_shamir_verifier(PROOF_BYTES)));
    }

    let last_codeword = match proof_stream.pull() {
      Some(StarkProofStreamEnum::Codeword(c)) => c,
      _ => {
        return Err(format!("expected {}th item in proof stream to be a codeword", num_rounds-1));
      }
    };

    if &MerkleRoot::commit(&last_codeword) != roots.last().unwrap() {
      return Err("last codeword is not well formed".to_string());
    }

    // check if it is low degree
    let degree = (last_codeword.len() / self.expansion_factor) - 1;
    let mut last_omega = omega;
    let mut last_offset = offset;
    for _ in 0..num_rounds-1 {
      last_omega = last_omega ^ 2_usize;
      last_offset = last_offset ^ 2_usize;
    }

    if last_omega.inverse() != (last_omega ^ (last_codeword.len() - 1)) {
      return Err("omega does not have the right order".to_string());
    }

    let last_domain = (0..last_codeword.len())
      .map(|i| last_offset * (last_omega ^ i))
      .collect::<Vec<_>>();
    let poly = Polynomial::interpolate_domain(&last_domain, &last_codeword);

    if poly.evaluate_domain(&last_domain) != last_codeword {
      return Err("re-evaluated codeword does not match original".to_string());
    }

    match poly.degree() {
      None => {
        return Err("Received none instead of polynomial degree".to_string());
      },
      Some(poly_degree) => {
        if poly_degree > degree {
          return Err(format!(
            "last codeword does not correspond to polynomial of low enough degree (it is {} but should be {})",
            poly_degree,
            degree,
          ));
        }
      }
    }

    let top_level_indices = self.sample_indices(
      proof_stream.fiat_shamir_verifier(PROOF_BYTES),
      self.domain_length >> 1,
      self.domain_length >> (num_rounds - 1),
      self.num_colinearity_tests,
    );

    for r in 0..num_rounds-1 {
      let indices_c = top_level_indices
        .iter()
        .map(|i| {
          i % (self.domain_length >> (r + 1))
        })
        .collect::<Vec<_>>();

      let indices_a = indices_c.clone();
      let indices_b = indices_a
        .iter()
        .map(|i| {
          i + (self.domain_length >> (r + 1))
        })
        .collect::<Vec<_>>();

      let mut aa = vec![];
      let mut bb = vec![];
      let mut cc = vec![];
      for s in 0..self.num_colinearity_tests {
        let (ay, by, cy) = match proof_stream.pull() {
          Some(StarkProofStreamEnum::Leafs(l)) => l,
          _ => {
            return Err(format!(
              "expected {}th item in proof stream to be a leaf",
              num_rounds + r * self.num_colinearity_tests * 4 + s,
            ));
          }
        };

        aa.push(ay);
        bb.push(by);
        cc.push(cy);

        if r == 0 {
          polynomial_values.push((*indices_a.get(s).unwrap(), ay));
          polynomial_values.push((*indices_b.get(s).unwrap(), by));
        }

        // colinearity check
        let ax = offset * (omega ^ *indices_a.get(s).unwrap());
        let bx = offset * (omega ^ *indices_b.get(s).unwrap());
        let cx = alphas[r];
        if !Polynomial::test_colinearity(vec![
          (ax, ay),
          (bx, by),
          (cx, cy),
        ]) {
          return Err("colinearity check failure".to_string());
        }
      }

      for i in 0..self.num_colinearity_tests {
        let nth = num_rounds + r * self.num_colinearity_tests * 4 + self.num_colinearity_tests + 3 * i;

        let path = match proof_stream.pull() {
          Some(StarkProofStreamEnum::Path(p)) => p,
          _ => {
            return Err(format!("expected {}th item in proof stream to be a path", nth + 1));
          }
        };
        if !MerkleRoot::verify(roots.get(r).unwrap(), *indices_a.get(i).unwrap(), &path, aa.get(i).unwrap()) {
          return Err("Merkle auth path verification failed for aa".to_string());
        }

        let path = match proof_stream.pull() {
          Some(StarkProofStreamEnum::Path(p)) => p,
          _ => {
            return Err(format!("expected {}th item in proof stream to be a path", nth + 2));
          }
        };
        if !MerkleRoot::verify(roots.get(r).unwrap(), *indices_b.get(i).unwrap(), &path, bb.get(i).unwrap()) {
          return Err("Merkle auth path verification failed for bb".to_string());
        }

        let path = match proof_stream.pull() {
          Some(StarkProofStreamEnum::Path(p)) => p,
          _ => {
            return Err(format!("expected {}th item in proof stream to be a path", nth + 3));
          }
        };
        if !MerkleRoot::verify(roots.get(r + 1).unwrap(), *indices_c.get(i).unwrap(), &path, cc.get(i).unwrap()) {
          return Err("Merkle auth path verification failed for cc".to_string());
        }
      }

      omega = omega ^ 2_usize;
      offset = offset ^ 2_usize;
    }

    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use crate::field::field::FIELD_PRIME;
  use super::*;

  #[test]
  fn sample_indices () {
    let n = 256;

    let field = Field::new(FIELD_PRIME);
    let fri = FRI::new(
      field.generator(),
      field.primitive_nth_root(n as u128),
      n,
      4,
      17,
    );

    let sample = fri.sample_indices(
      "d4b6e8af1114859c1c24b6496a3aef2f55a21105bc103af7e12dc3b2c101fe66".into(),
      128,
      128,
      17,
    );
    assert_eq!(
      sample,
      vec![40, 121, 5, 113, 97, 68, 126, 88, 26, 82, 81, 91, 93, 125, 10, 57, 48],
    )
  }

  #[test]
  fn verify () {
    // constants
    let field = Field::new(FIELD_PRIME);
    let degree = 63; // power of 2
    let expansion_factor = 4; // power of 2
    let num_colinearity_tests = 17;

    // calculated
    let codeword_initial_length = (degree + 1) * expansion_factor;
    let mut codeword_log_length = 0;
    let mut codeword_length = codeword_initial_length;
    while codeword_length > 1 {
      codeword_length /= 2;
      codeword_log_length += 1;
    }

    assert_eq!(1 << codeword_log_length, codeword_initial_length, "log calculated incorrectly");

    let omega = field.primitive_nth_root(codeword_initial_length as u128);
    let generator = field.generator();

    let fri = FRI::new(
      generator,
      omega,
      codeword_initial_length,
      expansion_factor,
      num_colinearity_tests,
    );

    let polynomial = Polynomial::new(
      (0..(degree as u128 + 1))
        .map(|i| FieldElement {
          field: &field,
          value: i,
        })
        .collect()
    );
    let domain = (0..codeword_initial_length)
      .map(|i| omega ^ i)
      .collect::<Vec<_>>();

    let codeword = polynomial.evaluate_domain(&domain);

    let mut proof_stream = StarkProofStream::new();

    println!("proving...");
    // todo: return proofstream
    fri.prove(&codeword, &mut proof_stream);
    println!("done");

    println!("{}", proof_stream.serialize());

    let mut points = vec![];
    println!("verifying...");
    // todo: consume proof_stream and return points
    assert_eq!(fri.verify(&mut proof_stream, &mut points), Ok(()), "proof should be valid");

    points
        .iter()
        .for_each(|(x, y)| {
          assert_eq!(&polynomial.evaluate(&(omega ^ *x)), y, "polynomial evaluates to wrong value");
        });

    println!("testing invalid codeword");

    // disturb then test for failure
    let mut codeword = codeword;
    for i in 0..degree/3 {
      codeword
          .get_mut(i)
          .map(|i| {
            *i = field.zero()
          });
    }
    let mut proof_stream = StarkProofStream::new();

    fri.prove(&codeword, &mut proof_stream);

    let mut points = vec![];
    assert_ne!(fri.verify(&mut proof_stream, &mut points), Ok(()), "proof should fail, but is accepted ");

    println!("success! \\o/");
  }
}
