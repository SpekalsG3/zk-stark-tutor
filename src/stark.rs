use std::fmt::{Display, Formatter};
use rand::{RngCore, thread_rng};
use serde::Serialize;
use crate::field::field::Field;
use crate::field::field_element::FieldElement;
use crate::field::polynomial::Polynomial;
use crate::fri::FRI;
use crate::m_polynomial::MPolynomial;
use crate::merkle_root::MerkleRoot;
use crate::proof_stream::PROOF_BYTES;
use crate::utils::bit_iter::BitIter;
use crate::utils::bytes::Bytes;

#[derive(Debug, Clone, Serialize)]
pub enum ProofStreamEnum<'a> {
    Root(Bytes),
    Codeword(Vec<FieldElement<'a>>),
    Path(Vec<Bytes>),
    Leafs((FieldElement<'a>, FieldElement<'a>, FieldElement<'a>)),
    Value(FieldElement<'a>),
}

impl Display for ProofStreamEnum<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            ProofStreamEnum::Root(r) => {
                write!(f, "\"{}\"", r.to_hex())
            }
            ProofStreamEnum::Codeword(c) => {
                let mut iter = c.iter();
                if let Some(o) = iter.next() {
                    write!(f, "[{}", o.value)?;
                }
                for o in iter {
                    write!(f, ",{}", o.value)?;
                }
                write!(f, "]")
            }
            ProofStreamEnum::Path(p) => {
                let mut iter = p.iter();
                if let Some(o) = iter.next() {
                    write!(f, "[{}", o.to_hex())?;
                }
                for o in iter {
                    write!(f, ",{}", o.to_hex())?;
                }
                write!(f, "]")
            }
            ProofStreamEnum::Leafs(l) => {
                write!(f,
                       "({},{},{})",
                       l.0.value,
                       l.1.value,
                       l.2.value,
                )
            }
            ProofStreamEnum::Value(fe) => {
                write!(f, "{}", fe.value)
            }
        }
    }
}

pub type StarkProofStream<'a> = crate::proof_stream::ProofStream<ProofStreamEnum<'a>>;

pub struct Stark<'a> {
    field: &'a Field,
    expansion_factor: usize,
    num_collinearity_checks: usize,
    security_level: usize,
    num_registers: usize,
    original_trace_length: usize,
    transition_constraints_degree: usize,
    num_randomizers: usize,
    omicron: FieldElement<'a>,
    omicron_domain: Vec<FieldElement<'a>>,
    fri: FRI<'a>,
}

impl<'a> Stark<'a> {
    pub fn new (
        field: &'a Field,
        expansion_factor: usize,
        num_collinearity_checks: usize,
        security_level: usize,
        num_registers: usize,
        num_cycles: usize,
        transition_constraints_degree: usize,
    ) -> Self {
        assert!(BitIter::from(field.order).count() >= security_level, "field order has to be at least {} bits", security_level);
        assert_eq!(expansion_factor & (expansion_factor - 1), 0, "expansion_factor must be a power of 2");
        assert!(expansion_factor >= 4, "expansion_factor must be at least 4");
        assert!(num_collinearity_checks * 2 >= security_level, "number of collinearity checks must be at least half of {}", security_level);

        let num_randomizers = 4 * num_collinearity_checks;
        let randomized_trace_length = num_cycles + num_randomizers;
        let omicron_domain_length = 1 << BitIter::from(randomized_trace_length * transition_constraints_degree).count();
        let fri_domain_length = omicron_domain_length * expansion_factor;

        let generator = field.generator();
        let omega = field.primitive_nth_root(fri_domain_length as u128);
        let omicron = field.primitive_nth_root(omicron_domain_length as u128);
        let omicron_domain = (0..omicron_domain_length)
            .map(|i| omicron ^ i as u128)
            .collect::<Vec<_>>();

        Self {
            field,
            expansion_factor,
            num_collinearity_checks,
            security_level,
            num_registers,
            original_trace_length: num_cycles,
            transition_constraints_degree,
            num_randomizers,
            omicron,
            omicron_domain,
            fri: FRI::new(
                generator,
                omega,
                fri_domain_length,
                expansion_factor,
                num_collinearity_checks
            ),
        }
    }

    pub fn transition_degree_bounds(
        &self,
        transition_constraints: &[MPolynomial],
    ) -> Vec<u128> {
        let mut points_degree = vec![1; 1];
        points_degree.resize(
            self.num_registers * 2 + 1,
            (self.original_trace_length + self.num_randomizers - 1) as u128
        );

        let mut res = vec![];
        for a in transition_constraints {
            if a.dictionary.len() == 0 {
                panic!("cannot calculate max on empty vec a");
            }

            let mut max = 0;
            for (k, _) in &a.dictionary {
                let mut sum = 0;

                let mut iter_points = points_degree.iter();
                let mut iter_k = k.iter();

                loop {
                    let r = match iter_points.next() {
                        None => break,
                        Some(r) => r,
                    };
                    let l = match iter_k.next() {
                        None => break,
                        Some(l) => l,
                    };
                    sum += r * l;
                }

                if sum > max {
                    max = sum
                }
            }

            res.push(max);
        }

        res
    }

    pub fn transition_quotient_degree_bounds(
        &self,
        transition_constraints: &[MPolynomial],
    ) -> Vec<u128> {
        self.transition_degree_bounds(transition_constraints)
            .iter()
            .map(|d| d - (self.original_trace_length - 1) as u128)
            .collect()
    }

    pub fn max_degree(
        &self,
        transition_constraints: &[MPolynomial],
    ) -> u128 {
        if transition_constraints.len() == 0 {
            panic!("Cannot calculate max for empty transition_constraints vector")
        }

        let md = self.transition_degree_bounds(transition_constraints)
            .into_iter()
            .max()
            .unwrap();

        (1 << BitIter::from(md).count()) - 1
    }

    pub fn transition_zerofier(&self) -> Polynomial {
        let domain = &self.omicron_domain[0..self.original_trace_length-1];
        Polynomial::zerofier_domain(domain)
    }

    pub fn boundary_zerofiers(
        &self,
        boundary: &[(usize, usize, FieldElement)],
    ) -> Vec<Polynomial> {
        (0..self.num_registers)
            .map(|s| {
                let domain = boundary
                    .iter()
                    .filter(|(_, r,_ )| r == &s)
                    .map(|(c, _, _)| self.omicron ^ *c)
                    .collect::<Vec<_>>();
                Polynomial::zerofier_domain(&domain)
            })
            .collect()
    }

    pub fn boundary_interpolants<'m>(
        &'m self,
        boundary: &[(usize, usize, FieldElement<'m>)],
    ) -> Vec<Polynomial> {
        (0..self.num_registers)
            .map(|s| {
                let mut domain = Vec::with_capacity(boundary.len());
                let mut values = Vec::with_capacity(boundary.len());

                for (c, r, v) in boundary {
                    if r != &s {
                        continue
                    }

                    domain.push(self.omicron ^ *c);
                    values.push(*v);
                }

                Polynomial::interpolate_domain(&domain, &values)
            })
            .collect()
    }

    pub fn boundary_quotient_degree_bounds(
        &self,
        randomized_trace_length: usize,
        boundary: &[(usize, usize, FieldElement)],
    ) -> Vec<usize> {
        let randomized_trace_degree = randomized_trace_length - 1;
        self.boundary_interpolants(boundary)
            .into_iter()
            .map(|bz| randomized_trace_degree - bz.degree().expect("Couldnt get degree of boundary zerofier"))
            .collect()
    }

    pub fn sample_weights(
        &self,
        number: usize,
        randomness: &Bytes,
    ) -> Vec<FieldElement> {
        (0..number)
            .map(|i| {
                self.field.sample(
                    &( Bytes::new(vec![0; i]) + randomness.clone() )
                )
            })
            .collect()
    }

    pub fn prove<'m>(
        &'m self,
        trace: Vec<Vec<FieldElement<'m>>>,
        transition_constraints: &[MPolynomial<'m>],
        boundary: &[(usize, usize, FieldElement<'m>)],
        proof_stream: Option<StarkProofStream>, // todo: Option here is questionable
    ) -> String {
        let mut thread_rng = thread_rng();
        let mut proof_stream = proof_stream.unwrap_or(StarkProofStream::new());

        // concatenate randomizers - induces zero-knowledge
        let trace = {
            let mut trace = trace;
            let trace_len = trace.len();

            trace.resize(trace_len + self.num_randomizers, Vec::with_capacity(self.num_registers));
            for i in 0..self.num_randomizers {
                trace[trace_len + i]
                    .fill_with(|| {
                        let mut bytes = vec![0; 17]; // todo: shouldn't 17 be self.num_colinearity_tests?
                        thread_rng.fill_bytes(&mut bytes);
                        self.field.sample(&Bytes::new(bytes))
                    })
            };

            trace
        };

        // interpolate
        // AIR = arithmetic constraint system = execution trace = list of trace polynomials - proof is valid iff AIR has a satisfying solution
        let trace_polynomials = {
            let trace_domain = (0..trace.len())
                .map(|i| self.omicron ^ i)
                .collect::<Vec<_>>();

            (0..self.num_registers)
                .map(|s| {
                    let single_trace = trace
                        .iter()
                        .map(|el| *el.get(s).expect("One of trace elements' length is less then self.num_registers"))
                        .collect::<Vec<_>>();

                    Polynomial::interpolate_domain(&trace_domain, &single_trace)
                })
                .collect::<Vec<_>>()
        };

        // subtract boundary interpolants and divide out boundary zerofiers
        // both quotients are required for the verifier to check their bounded degree (what FRI solves)
        // additionally verifier checks link between both quotients
        let boundary_quotients = {
            let boundary_interpolants = self.boundary_interpolants(boundary);
            let boundary_zerofiers = self.boundary_zerofiers(boundary);

            (0..self.num_registers)
                .map(|s| {
                    let interpolant = boundary_interpolants
                        .get(s)
                        .expect("Number of boundary_interpolants is less then self.num_registers")
                        .clone();
                    let zerofier = boundary_zerofiers
                        .get(s)
                        .expect("Number of boundary_zerofiers is less then self.num_registers")
                        .clone();
                    // SAFETY: [s] safe because trace_polynomials also size of self.num_registers
                    let trace_polynomial = trace_polynomials[s].clone();
                    let boundary_polynomial = trace_polynomial - interpolant;

                    boundary_polynomial / zerofier
                })
                .collect::<Vec<_>>()
        };

        // commit to boundary_quotients
        // trace polynomials relate to both quotients (but not equivalent to transition quotient, which doesn't matter)
        // and thus we don't need to commit to both quotients
        // and we don't need to commit to trace polynomials because if verifier knows a value in a given
        // point of one quotient, he can compute a matching value of another using only public data
        let fri_domain = self.fri.eval_domain();
        let boundary_quotient_codewords = {
            (0..self.num_registers)
                .map(|s| {
                    // SAFETY: [s] safe because trace_polynomials also size of self.num_registers
                    let codeword = boundary_quotients[s].clone().evaluate_domain(&fri_domain);

                    let root = MerkleRoot::commit(&codeword);
                    proof_stream.push(ProofStreamEnum::Root(root));

                    codeword
                })
                .collect::<Vec<_>>()
        };

        let transition_quotients = {
            let mut point = Vec::with_capacity(1 + trace_polynomials.len() * 2);
            point.push(Polynomial::new(vec![
                self.field.zero(),
                self.field.one(),
            ]));
            point.extend(trace_polynomials.clone());
            point.extend(
                trace_polynomials
                    .iter()
                    .map(|tp| tp.scale(self.omicron))
            );

            transition_constraints
                .iter()
                .map(|a| {
                    // symbolically evaluate transition constraints to receive transition_polynomial
                    let transition_polynomial = a.evaluate_symbolic(&point);

                    // divide out zerofier
                    transition_polynomial / self.transition_zerofier()
                })
                .collect::<Vec<_>>()
        };

        // commit to randomizer polynomial
        let randomizer_polynomial = Polynomial::new(
            (0..self.max_degree(&transition_constraints) + 1)
                .map(|_| {
                    let mut bytes = vec![0; 17]; // todo: shouldn't 17 be self.num_colinearity_tests?
                    thread_rng.fill_bytes(&mut bytes);
                    self.field.sample(&Bytes::new(bytes))
                })
                .collect()
        );
        let randomizer_codeword = randomizer_polynomial.evaluate_domain(&fri_domain);
        {
            let randomizer_root = MerkleRoot::commit(&randomizer_codeword);
            proof_stream.push(ProofStreamEnum::Root(randomizer_root));
        }

        let weights = self.sample_weights(
            1 + 2 * transition_quotients.len() + 2 * boundary_quotients.len(),
            &proof_stream.fiat_shamir_prover(PROOF_BYTES),
        );
        assert_eq!(
            transition_quotients.iter().map(|tq| tq.degree().expect("Failed to get degree of transition quotient") as u128).collect::<Vec<_>>(),
            self.transition_quotient_degree_bounds(&transition_constraints),
            "transition quotient degrees do not match with expectation",
        );

        let terms = {
            let x = Polynomial::new(vec![self.field.zero(), self.field.one()]);

            let mut res = Vec::with_capacity(weights.len());
            res.push(randomizer_polynomial);

            let transition_constraints_degree = self.max_degree(&transition_constraints);

            let transition_quotient_degree_bounds = self.transition_quotient_degree_bounds(&transition_constraints);
            for (i, tq) in transition_quotients.iter().enumerate() {
                res.push(tq.clone());

                // SAFETY: [i] safe because it's the same size as `transition_quotients`
                let shift = transition_constraints_degree - transition_quotient_degree_bounds[i];
                res.push((x.clone() ^ shift) * tq.clone());
            }

            let boundary_quotient_degree_bounds = self.boundary_quotient_degree_bounds(trace.len(), &boundary);
            for (i, bq) in boundary_quotients.iter().enumerate() {
                res.push(bq.clone());

                // SAFETY: [i] safe because it's the same size as `boundary_quotients`
                let shift = transition_constraints_degree - boundary_quotient_degree_bounds[i] as u128;
                res.push((x.clone() ^ shift) * bq.clone());
            }

            res
        };

        let indices = {
            // take weighted sum
            let combination = terms
                .into_iter()
                .enumerate()
                .map(|(i, term)| {
                    // SAFETY: [i] safe because weights and terms have the same size
                    Polynomial::new(vec![ weights[i] ]) * term
                })
                .reduce(|a, b| a + b)
                .unwrap(); // SAFETY: terms length is at least 1

            // compute matching codeword
            let combined_codeword = combination.evaluate_domain(&fri_domain);

            // prove low degree of combination polynomial
            let mut indices = self.fri.prove(&combined_codeword, &mut proof_stream);
            indices.sort();
            indices
        };
        let duplicated_indices = {
            let mut v = Vec::with_capacity(2 * indices.len());
            v.extend(indices.clone());
            v.extend(
                indices
                    .iter()
                    .map(|i| {
                        (i + self.expansion_factor) % fri_domain.len()
                    })
            );
            v
        };

        // open indicated positions in the boundary quotient codewords
        for bqc in boundary_quotient_codewords {
            for i in &duplicated_indices {
                // SAFETY: [i] is safe because bqc is size of fri_domain_length and it's much bigger than size of duplicated_indices (2*num_collinearity_checks)
                proof_stream.push(ProofStreamEnum::Value(bqc[*i]));
                let path = MerkleRoot::open(*i, &bqc);
                proof_stream.push(ProofStreamEnum::Path(path));
            }
        }
        // and in the randomizer
        for i in indices {
            // SAFETY: [i] safe because fri_domain_length is much bigger then size of indices (num_collinearity_checks)
            proof_stream.push(ProofStreamEnum::Value(randomizer_codeword[i]));
            let path = MerkleRoot::open(i, &randomizer_codeword);
            proof_stream.push(ProofStreamEnum::Path(path));
        }

        proof_stream.serialize()
    }

    pub fn verify<'m> (
        &self,
        proof: String,
        transition_constraints: &[MPolynomial<'m>],
        boundary: &[(usize, usize, FieldElement<'m>)],
        proof_stream: Option<StarkProofStream>, // todo: Option here is questionable
    ) -> bool {
        let mut proof_stream = proof_stream.unwrap_or(StarkProofStream::new());
        unimplemented!()
    }
}
