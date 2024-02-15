use std::io::Read;
use crate::crypto::shake256::shake256;
use crate::utils::bytes::Bytes;
use crate::utils::stringify::Stringify;

#[derive(Debug)]
pub struct ProofStream<T> {
  objects: Vec<T>,
  read_index: usize,
}

impl<T: PartialEq> PartialEq for ProofStream<T> {
  fn eq(&self, other: &Self) -> bool {
    self.objects.iter().eq(other.objects.iter())
  }
}

impl<T> ProofStream<T>
  where
    for<'a> &'a[T]: Stringify,
{
  pub fn serialize (&self) -> Bytes {
    self.objects.as_slice().stringify().as_bytes().into()
  }

  // get challenge
  pub fn fiat_shamir_prover (&self, num_bytes: usize) -> Bytes {
    let str = self.serialize();

    shake256(str, num_bytes)
  }

  // reproduce challenge
  pub fn fiat_shamir_verifier (&self, num_bytes: usize) -> Bytes {
    let slice = &self.objects[0..self.read_index];
    let str = slice.stringify();

    shake256(str.as_bytes().into(), num_bytes)
  }
}

impl<T> ProofStream<T>
  where
    T: Clone,
    for<'a> &'a[T]: Stringify,
{
  pub fn new () -> Self {
    ProofStream {
      objects: vec![],
      read_index: 0,
    }
  }

  pub fn from (objects: Vec<T>) -> Self {
    Self {
      objects,
      read_index: 0,
    }
  }

  pub fn push (&mut self, obj: T) {
    self.objects.push(obj);
  }

  pub fn pull (&mut self) -> Option<T> {
    assert!(self.read_index < self.objects.len(), "Cannot pull, queue is empty");

    let obj = match self.objects.get(self.read_index) {
      Some(o) => Some(o.clone()),
      None => None
    };
    self.read_index += 1;

    obj
  }
}

#[cfg(test)]
mod tests {
  use std::collections::HashMap;
  use crate::proof_stream::ProofStream;
  use crate::utils::stringify::Stringify;

  #[test]
  fn order () {
    #[derive(Debug, Clone)]
    enum SomeObjects {
      Map(HashMap<String, usize>),
      Vec(Vec<usize>),
      Str(String),
    }
    impl PartialEq for SomeObjects {
      fn eq (&self, other: &Self) -> bool {
        match self {
          SomeObjects::Map(m_1) => {
            match other {
              SomeObjects::Map(m_2) => m_1.eq(m_2),
              _ => false,
            }
          },
          SomeObjects::Vec(v_1) => {
            match other {
              SomeObjects::Vec(v_2) => v_1.eq(v_2),
              _ => false,
            }
          },
          SomeObjects::Str(s_1) => {
            match other {
              SomeObjects::Str(s_2) => s_1.eq(s_2),
              _ => false,
            }
          },
        }
      }
    }
    impl Stringify for &[SomeObjects] {
      fn stringify<'m>(&'m self) -> String {
        format!("{:?}", self) // in test we don't care about optimized size
      }
    }

    let mut proof_stream = ProofStream::<SomeObjects>::new();

    assert_eq!(proof_stream.fiat_shamir_prover(64).to_hex(), "ec784925b52067bce01fd820f554a34a3f8522b337f82e00ea03d3fa2b207ef9c2c1b9ed900cf2bbfcd19a232a94c6121e041615305c4155d46d52f58a8cff1c");

    proof_stream.push(SomeObjects::Str(String::from("Hello, World!")));
    proof_stream.push(SomeObjects::Vec(vec![0, 1, 5, 234]));
    let mut map = HashMap::new();
    map.insert(String::from("something"), 123);
    proof_stream.push(SomeObjects::Map(map.clone()));

    assert_eq!(proof_stream.fiat_shamir_prover(4).to_hex(), "78b0db5c");
    assert_eq!(proof_stream.fiat_shamir_prover(64).to_hex(), "78b0db5cfd13c78498fd0951a9fd609f2521fd02d850cc561eced844bb0c338588358abcc0d98d76c6779cb388514f4bc19e2c0125b143abee166cb98c38a831");
    assert_eq!(proof_stream.fiat_shamir_verifier(64).to_hex(), "ec784925b52067bce01fd820f554a34a3f8522b337f82e00ea03d3fa2b207ef9c2c1b9ed900cf2bbfcd19a232a94c6121e041615305c4155d46d52f58a8cff1c");

    assert_eq!(proof_stream.pull(), Some(SomeObjects::Str(String::from("Hello, World!"))));
    assert_eq!(proof_stream.pull(), Some(SomeObjects::Vec(vec![0, 1, 5, 234])));
    assert_eq!(proof_stream.pull(), Some(SomeObjects::Map(map)));

    assert_eq!(proof_stream.fiat_shamir_verifier(64).to_hex(), "78b0db5cfd13c78498fd0951a9fd609f2521fd02d850cc561eced844bb0c338588358abcc0d98d76c6779cb388514f4bc19e2c0125b143abee166cb98c38a831");
  }
}
