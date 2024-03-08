use std::collections::HashMap;
use std::io::Read;
use rand::{RngCore, thread_rng};
use crate::crypto::shake256::PROOF_BYTES;
use crate::fft::ntt_arithmetics::{fast_coset_divide, fast_coset_evaluate, fast_interpolate_domain, fast_multiply, fast_zerofier};
use crate::field::field::Field;
use crate::field::field_element::FieldElement;
use crate::field::polynomial::Polynomial;
use crate::fri::FRI;
use crate::m_polynomial::MPolynomial;
use crate::merkle_root::MerkleRoot;
use crate::proof_stream::{IndependentProofStream, ProofStream};
use crate::stark::proof_stream_enum::StarkProofStreamEnum;
use crate::utils::bit_iter::BitIter;
use crate::utils::bytes::Bytes;

pub struct Stark<'a> {
    pub(crate) field: &'a Field,
    expansion_factor: usize,
    num_registers: usize,
    original_trace_length: usize,
    num_randomizers: usize,
    pub(crate) omicron: FieldElement<'a>,
    omicron_domain: Vec<FieldElement<'a>>,
    pub(crate) omicron_domain_length: u128,
    fri: FRI<'a>,
}

impl<'a> Stark<'a> {
    pub fn deser_independent_proof_stream(&self, bytes: Bytes) -> IndependentProofStream<StarkProofStreamEnum> {
        let mut b = bytes.bytes();

        {
            let mut order = [0; 128 / 8];
            b.read_exact(&mut order).unwrap();
            let order = u128::from_be_bytes(order);
            if order != 0 {
                if order != self.field.order {
                    panic!("serialized field differs from Stark's field");
                }
            }
        };

        let mut items = vec![];
        loop {
            let mut code = [0; 1];
            if let Err(_) = b.read_exact(&mut code) {
                break
            };

            let el = {
                let mut size = [0; 64 / 8];
                b.read_exact(&mut size).unwrap();
                let size = usize::from_be_bytes(size);

                let mut el = vec![0; size];
                b.read_exact(&mut el).unwrap();

                Bytes::new(el)
            };

            let item = StarkProofStreamEnum::from_bytes(code[0], el, self.field);
            items.push(item)
        }

        IndependentProofStream::from(items)
    }
}

impl<'a> Stark<'a> {
    pub fn new (
        field: &'a Field,
        expansion_factor: usize,
        num_collinearity_checks: usize,
        security_level: usize,
        num_registers: usize,
        num_cycles: usize,
        transition_constraints_degree: usize, // 2
    ) -> Self {
        assert!(BitIter::from(field.order).count() >= security_level, "field order has to be at least {} bits", security_level);
        assert_eq!(expansion_factor & (expansion_factor - 1), 0, "expansion_factor must be a power of 2");
        assert!(expansion_factor >= 4, "expansion_factor must be at least 4");
        assert!(num_collinearity_checks * 2 >= security_level, "number of collinearity checks must be at least half of {}", security_level);

        let num_randomizers = 4 * num_collinearity_checks;
        let randomized_trace_length = num_cycles + num_randomizers;
        let omicron_domain_length = 1_u128 << BitIter::from(randomized_trace_length * transition_constraints_degree).count();
        let fri_domain_length = omicron_domain_length * (expansion_factor as u128);

        let generator = field.generator();
        let omega = field.primitive_nth_root(fri_domain_length);
        let omicron = field.primitive_nth_root(omicron_domain_length);
        let omicron_domain = (0..omicron_domain_length)
            .map(|i| omicron ^ i)
            .collect::<Vec<_>>();

        Self {
            field,
            expansion_factor,
            num_registers,
            original_trace_length: num_cycles,
            num_randomizers,
            omicron,
            omicron_domain,
            omicron_domain_length,
            fri: FRI::new(
                generator,
                omega,
                fri_domain_length as usize,
                expansion_factor,
                num_collinearity_checks
            ),
        }
    }

    fn transition_degree_bounds(
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

    fn transition_quotient_degree_bounds(
        &self,
        transition_constraints: &[MPolynomial],
    ) -> Vec<u128> {
        self.transition_degree_bounds(transition_constraints)
            .iter()
            .map(|d| d - (self.original_trace_length - 1) as u128)
            .collect()
    }

    fn max_degree(
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

    fn transition_zerofier(&self) -> Polynomial {
        let domain = &self.omicron_domain[0..self.original_trace_length-1];
        // Polynomial::zerofier_domain(domain)
        fast_zerofier(
            self.omicron,
            self.omicron_domain_length,
            &domain,
        )
    }

    fn boundary_zerofiers(
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
                // Polynomial::zerofier_domain(&domain)
                fast_zerofier(
                    self.omicron,
                    self.omicron_domain_length,
                    &domain,
                )
            })
            .collect()
    }

    fn boundary_interpolants<'m>(
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

                // Polynomial::interpolate_domain(&domain, &values)
                fast_interpolate_domain(
                    self.omicron,
                    self.omicron_domain_length,
                    &domain,
                    &values
                )
            })
            .collect()
    }

    fn boundary_quotient_degree_bounds(
        &self,
        randomized_trace_length: usize,
        boundary: &[(usize, usize, FieldElement)],
    ) -> Vec<usize> {
        let randomized_trace_degree = randomized_trace_length - 1;
        self.boundary_zerofiers(boundary)
            .into_iter()
            .map(|bz| {
                randomized_trace_degree - bz.degree().expect("Couldnt get degree of boundary zerofier")
            })
            .collect()
    }

    fn sample_weights(
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

    pub fn prove<'m, PS: ProofStream<StarkProofStreamEnum<'m>>>(
        &'m self,
        trace: Vec<Vec<FieldElement<'m>>>,
        transition_constraints: &[MPolynomial<'m>],
        boundary: &[(usize, usize, FieldElement<'m>)],
        mut proof_stream: PS,
    ) -> Result<Bytes, String> {
        let mut thread_rng = thread_rng();

        // concatenate randomizers - induces zero-knowledge
        let trace = {
            let mut trace = trace;
            let trace_len = trace.len();

            trace.resize_with(trace_len + self.num_randomizers, || {
                (0..self.num_registers)
                    .map(|_| {
                        let mut bytes = vec![0; 17]; // todo: shouldn't 17 be self.num_colinearity_tests?
                        thread_rng.fill_bytes(&mut bytes);
                        self.field.sample(&Bytes::new(bytes))
                    })
                    .collect()
            });

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

                    // Polynomial::interpolate_domain(&trace_domain, &single_trace)
                    fast_interpolate_domain(
                        self.omicron,
                        self.omicron_domain_length,
                        &trace_domain,
                        &single_trace,
                    )
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

                    // boundary_polynomial.div(zerofier)
                    let q = fast_coset_divide(
                        self.omicron,
                        self.omicron_domain_length,
                        self.field.generator(),
                        boundary_polynomial,
                        zerofier,
                    );
                    Ok::<_, String>(q)
                })
                .collect::<Result<Vec<_>, _>>()?
        };

        // commit to boundary_quotients
        // trace polynomials relate to both quotients (but not equivalent to transition quotient, which doesn't matter)
        // and thus we don't need to commit to both quotients
        // and we don't need to commit to trace polynomials because if verifier knows a value in a given
        // point of one quotient, he can compute a matching value of another using only public data
        let fri_domain = self.fri.evaluate_domain();
        let boundary_quotient_codewords = {
            (0..self.num_registers)
                .map(|s| {
                    // SAFETY: [s] safe because boundary_quotients also size of self.num_registers
                    // let codeword = boundary_quotients[s].evaluate_domain(&fri_domain);
                    let codeword = fast_coset_evaluate(
                        self.fri.omega,
                        self.fri.domain_length as u128,
                        self.field.generator(),
                        boundary_quotients[s].clone(),
                    );

                    let root = MerkleRoot::commit(&codeword);
                    proof_stream.push(StarkProofStreamEnum::Root(root));

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
                    .into_iter()
                    .map(|tp| tp.scale(self.omicron))
            );

            transition_constraints
                .iter()
                .map(|a| {
                    // symbolically evaluate transition constraints to receive transition_polynomial
                    let transition_polynomial = a.evaluate_symbolic(&point);

                    // divide out zerofier
                    // transition_polynomial.div(self.transition_zerofier())
                    let q = fast_coset_divide(
                        self.omicron,
                        self.omicron_domain_length,
                        self.field.generator(),
                        transition_polynomial,
                        self.transition_zerofier(),
                    );
                    Ok::<_, String>(q)
                })
                .collect::<Result<Vec<_>, _>>()?
        };

        // commit to randomizer polynomial
        let transition_constraints_degree = self.max_degree(&transition_constraints);

        let randomizer_polynomial = Polynomial::new(
            (0..(transition_constraints_degree + 1))
                .map(|_| {
                    let mut bytes = vec![0; 17]; // todo: shouldn't 17 be self.num_colinearity_tests?
                    thread_rng.fill_bytes(&mut bytes);
                    self.field.sample(&Bytes::new(bytes))
                })
                .collect()
        );

        // let randomizer_codeword = randomizer_polynomial.evaluate_domain(&fri_domain);
        let randomizer_codeword = fast_coset_evaluate(
            self.fri.omega,
            self.fri.domain_length as u128,
            self.field.generator(),
            randomizer_polynomial.clone(),
        );
        {
            let randomizer_root = MerkleRoot::commit(&randomizer_codeword);
            proof_stream.push(StarkProofStreamEnum::Root(randomizer_root));
        }

        let weights = self.sample_weights(
            1 + 2 * transition_quotients.len() + 2 * boundary_quotients.len(),
            &proof_stream.fiat_shamir_prover(PROOF_BYTES),
        );
        if
            transition_quotients
                .iter()
                .map(|tq| tq
                    .degree()
                    .map_or_else(
                        || Err("Failed to get degree of transition quotient"),
                        |d| Ok(d as u128),
                    )
                )
                .collect::<Result<Vec<_>, _>>()?
            != self.transition_quotient_degree_bounds(&transition_constraints) {
            return Err("transition quotient degrees do not match with expectation".to_string())
        }

        let terms = {
            let x = Polynomial::new(vec![self.field.zero(), self.field.one()]);

            let mut terms = Vec::with_capacity(weights.len());
            terms.push(randomizer_polynomial);

            let transition_quotient_degree_bounds = self.transition_quotient_degree_bounds(&transition_constraints);
            for (i, tq) in transition_quotients.iter().enumerate() {
                terms.push(tq.clone());

                // SAFETY: [i] safe because it's the same size as `transition_quotients`
                let shift = transition_constraints_degree - transition_quotient_degree_bounds[i];

                // let t = (x.clone() ^ shift) * tq.clone();
                let t = fast_multiply(self.omicron, self.omicron_domain_length, x.clone() ^ shift, tq.clone());

                terms.push(t);
            }

            let boundary_quotient_degree_bounds = self.boundary_quotient_degree_bounds(trace.len(), &boundary);
            for (i, bq) in boundary_quotients.into_iter().enumerate() {
                terms.push(bq.clone());

                // SAFETY: [i] safe because it's the same size as `boundary_quotients`
                let shift = transition_constraints_degree - boundary_quotient_degree_bounds[i] as u128;

                // let t = (x.clone() ^ shift) * bq.clone();
                let t = fast_multiply(self.omicron, self.omicron_domain_length, x.clone() ^ shift, bq);

                terms.push(t);
            }
            terms
        };

        let quadrupled_indices = {
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
            // let combined_codeword = combination.evaluate_domain(&fri_domain);
            let combined_codeword = fast_coset_evaluate(
                self.fri.omega,
                self.fri.domain_length as u128,
                self.field.generator(),
                combination,
            );

            // prove low degree of combination polynomial
            let indices = self.fri.prove(&combined_codeword, &mut proof_stream);

            let duplicated_indices = indices
                .clone()
                .into_iter()
                .chain(
                    indices
                        .iter()
                        .map(|i| {
                            (i + self.expansion_factor) % fri_domain.len()
                        })
                );
            let mut quadrupled_indices = duplicated_indices
                .clone()
                .chain(
                    duplicated_indices
                        .map(|i| (i + (fri_domain.len() / 2)) % fri_domain.len())
                )
                .collect::<Vec<_>>();
            quadrupled_indices.sort();
            quadrupled_indices
        };

        // open indicated positions in the boundary quotient codewords
        for bqc in boundary_quotient_codewords {
            for i in quadrupled_indices.clone() { // has to copy i anyhow so it's more efficient to clone whole vector
                // SAFETY: [i] is safe because bqc is size of fri_domain_length and it's much bigger than size of duplicated_indices (2*num_collinearity_checks)
                proof_stream.push(StarkProofStreamEnum::Value(bqc[i]));
                let path = MerkleRoot::open(i, &bqc);
                proof_stream.push(StarkProofStreamEnum::Path(path));
            }
        }
        // and in the randomizer
        for i in quadrupled_indices {
            // SAFETY: [i] safe because fri_domain_length is much bigger then size of indices (num_collinearity_checks)
            proof_stream.push(StarkProofStreamEnum::Value(randomizer_codeword[i]));
            let path = MerkleRoot::open(i, &randomizer_codeword);
            proof_stream.push(StarkProofStreamEnum::Path(path));
        }

        Ok(proof_stream.digest())
    }

    pub fn verify<'m, PS: ProofStream<StarkProofStreamEnum<'m>>> (
        &'m self,
        transition_constraints: &[MPolynomial<'m>],
        boundary: &[(usize, usize, FieldElement<'m>)],
        mut proof_stream: PS,
    ) -> Result<(), String> {
        // infer trace length from boundary conditions
        let original_trace_length = 1 + boundary
            .iter()
            .max_by(|(c1, _, _), (c2, _, _)| c1.cmp(c2))
            .expect("boundary is empty")
            .0;
        let randomized_trace_length = original_trace_length + self.num_randomizers;

        // get merkle roots of boundary quotient codewords
        let boundary_quotient_roots = (0..self.num_registers)
            .map(|_| proof_stream.pull().unwrap().expect_root())
            .collect::<Vec<_>>();

        // get merkle root of randomizer polynomial
        let randomizer_root = proof_stream.pull().unwrap().expect_root();

        let boundary_interpolants = self.boundary_interpolants(boundary);

        // get weights for non-linear combination
        let weights = self.sample_weights(
            1 + 2 * transition_constraints.len() + 2 * boundary_interpolants.len(),
            &proof_stream.fiat_shamir_verifier(PROOF_BYTES),
        );

        // verify low degree of combination polynomial
        let (
            indices,
            values,
        ) = {
            let mut points = vec![];
            if let Err(e) = self.fri.verify(
                &mut proof_stream,
                &mut points,
            ) {
                return Err(format!("FRI verification failed: {}", e));
            }
            points.sort_by(|a, b| a.0.cmp(&b.0));

            let points_len = points.len();
            points
                .into_iter().
                fold((
                    Vec::with_capacity(points_len),
                    Vec::with_capacity(points_len),
                ), |mut acc, p| {
                    acc.0.push(p.0);
                    acc.1.push(p.1);
                    acc
                })
        };
        let duplicated_indices = {
            let mut duplicated_indices = indices
                .clone()
                .into_iter()
                .chain(
                    indices
                        .iter()
                        .map(|i| {
                            (i + self.expansion_factor) % self.fri.domain_length
                        })
                )
                .collect::<Vec<_>>();
            duplicated_indices.sort();
            duplicated_indices
        };

        // read and verify leafs, which are elements of boundary quotient codewords
        let leafs = {
            boundary_quotient_roots
                .iter()
                .map(|bqr| {
                    duplicated_indices
                        .iter()
                        .map(|i| {
                            let leaf = proof_stream.pull().unwrap().expect_value();
                            let path = proof_stream.pull().unwrap().expect_path();
                            let accepted = MerkleRoot::verify(bqr, *i, &path, &leaf);

                            if !accepted {
                                return Err(format!("Boundary quotient root {} is not verified", i));
                            }

                            Ok((*i, leaf))
                        })
                        .collect::<Result<HashMap<_, _>, _>>()
                })
                .collect::<Result<Vec<_>, _>>()?
        };

        // read and verify randomizer leafs
        let randomizers = duplicated_indices
            .iter()
            .map(|i| {
                let leaf = proof_stream.pull().unwrap().expect_value();
                let path = proof_stream.pull().unwrap().expect_path();
                let accepted = MerkleRoot::verify(&randomizer_root, *i, &path, &leaf);

                if !accepted {
                    return Err(format!("Randomizer leaf {} not verified", i));
                }

                Ok((*i, leaf))
            })
            .collect::<Result<HashMap<_, _>, _>>()?;

        // verify leafs of combination polynomial
        let boundary_zerofiers = self.boundary_zerofiers(boundary);
        let boundary_interpolants = self.boundary_interpolants(boundary);
        indices
            .into_iter()
            .enumerate()
            .try_for_each(|(index_i, index_current)| {
                // get trace values by applying correction to boundary quotient values (which are leafs)
                let domain_index_current = self.fri.offset * (self.fri.omega ^ index_current);
                let index_next = (index_current + self.expansion_factor) % self.fri.domain_length;
                let domain_index_next = self.fri.offset * (self.fri.omega ^ index_next);

                let (
                    trace_current,
                    trace_next,
                ) = (0..self.num_registers)
                    .fold((
                        Vec::with_capacity(self.num_registers),
                        Vec::with_capacity(self.num_registers),
                    ), |mut acc, s| {
                        let zerofier = &boundary_zerofiers[s];
                        let interpolant = &boundary_interpolants[s];

                        acc.0.push(leafs[s][&index_current] * zerofier.evaluate(&domain_index_current) + interpolant.evaluate(&domain_index_current));
                        acc.1.push(leafs[s][&index_next]    * zerofier.evaluate(&domain_index_next)    + interpolant.evaluate(&domain_index_next));

                        acc
                    });

                let point = {
                    let mut v = Vec::with_capacity(1 + trace_current.len() + trace_next.len());
                    v.push(domain_index_current);
                    v.extend(trace_current.clone());
                    v.extend(trace_next.clone());
                    v
                };

                let transition_constraint_values = transition_constraints
                    .iter()
                    .map(|tc| tc.evaluate(&point))
                    .collect::<Vec<_>>();

                // compute non-linear combination
                let terms = {
                    let transition_constraints_degree = self.max_degree(transition_constraints);

                    let mut terms = Vec::with_capacity(1 + 2 * transition_constraint_values.len() + 2 * leafs.len());
                    terms.push(randomizers[&index_current]);

                    {
                        let transition_quotient_degree_bounds = self.transition_quotient_degree_bounds(transition_constraints);
                        // todo: how to name that?
                        let transition_quotient = self.transition_zerofier().evaluate(&domain_index_current);
                        transition_constraint_values
                            .into_iter()
                            .enumerate()
                            .for_each(|(s, tcv)| {
                                let quotient = tcv / transition_quotient;
                                terms.push(quotient);

                                let shift = transition_constraints_degree - transition_quotient_degree_bounds[s];
                                terms.push(quotient * (domain_index_current ^ shift));
                            });
                    }

                    {
                        let boundary_quotient_degree_bounds = self.boundary_quotient_degree_bounds(randomized_trace_length, boundary);
                        (0..self.num_registers)
                            .for_each(|s| {
                                let bqv = leafs[s][&index_current];
                                terms.push(bqv);

                                let shift = transition_constraints_degree - boundary_quotient_degree_bounds[s] as u128;
                                terms.push(bqv * (domain_index_current ^ shift));
                            });
                    }

                    terms
                };

                let combination = terms
                    .into_iter()
                    .enumerate()
                    .map(|(j, term)| term * weights[j])
                    .reduce(|a, b| a + b)
                    .unwrap(); // SAFETY: terms are at least length of 1

                // verify againts combination polynomial value
                if combination != values[index_i] {
                    return Err("Combination doesn't match with polynomial value".to_string())
                }

                Ok(())
            })
    }
}

#[cfg(test)]
mod tests {
    use rand::{RngCore, thread_rng};
    use crate::field::field::{Field, FIELD_PRIME};
    use crate::field::field_element::FieldElement;
    use crate::proof_stream::{IndependentProofStream, ProofStream};
    use crate::rescue_prime::rescue_prime::RescuePrime;
    use crate::stark::proof_stream_enum::StarkProofStreamEnum;
    use crate::stark::stark::Stark;
    use crate::utils::bytes::Bytes;

    #[test]
    fn deserialize_proof_stream () {
        let field = Field::new(FIELD_PRIME);
        let stark = Stark::new(
            &field,
            4, // expansion_factor,
            2, // num_collinearity_checks,
            2, // security_level,
            2, // num_registers,
            2, // num_cycles,
            2, // transition_constraints_degree,
        );
        let stream = IndependentProofStream::from(vec![
            StarkProofStreamEnum::Root(Bytes::new(vec![0x49,0x6e,0x20,0x74])),
            StarkProofStreamEnum::Codeword(vec![FieldElement::new(&field, 20), FieldElement::new(&field, 100)]),
            StarkProofStreamEnum::Path(vec![Bytes::new(vec![0x49,0x6e,0x20,0x74]), Bytes::new(vec![0x1,0x6b,0xfe,0x25])]),
            StarkProofStreamEnum::Leafs((FieldElement::new(&field, 1), FieldElement::new(&field, 5), FieldElement::new(&field, 10))),
            StarkProofStreamEnum::Value(FieldElement::new(&field, 2)),
        ]);

        let serialized = stream.digest();

        let deserialized = stark.deser_independent_proof_stream(serialized);
        assert_eq!(stream, deserialized);
    }

    #[test]
    fn stark () {
        let field = Field::new(FIELD_PRIME);
        let security_level = 2;
        let rp = RescuePrime::new(&field, 2, 1, security_level, 27);
        let stark = Stark::new(
            &field,
            4, // expansion_factor,
            2, // num_collinearity_checks,
            security_level,
            rp.m, // num_registers,
            rp.N + 1, // num_cycles,
            2, // transition_constraints_degree,
        );

        let mut output_element = field.sample(&Bytes::from("deadbeef"));
        let mut trace = vec![];
        let mut air = vec![];
        let mut boundary = vec![];

        for trial in 0..20 {
            let input_element = output_element.clone();
            println!("running trial #{} with input: {:?}", trial, input_element);
            output_element = rp.hash(input_element.clone());

            // prove honestly
            println!("honest proof generation ...");

            // prove
            trace = rp.trace(input_element);
            air = rp.transition_constraints(stark.omicron, stark.omicron_domain_length);
            boundary = rp.boundary_constraints(output_element);

            let proof = stark.prove(trace.clone(), &air, &boundary, IndependentProofStream::new());
            assert!(proof.is_ok(), "Failed to construct proof - {}", proof.unwrap_err());
            let proof = proof.unwrap();

            // verify
            let verdict = stark.verify(&air, &boundary, stark.deser_independent_proof_stream(proof.clone()));

            assert!(verdict.is_ok(), "valid stark proof fails to verify - {}", verdict.unwrap_err());
            println!("success \\o/");

            println!("verifying false claim ...");
            // verify false claim
            let output_element = output_element + field.one();
            let boundary = rp.boundary_constraints(output_element);
            let verdict = stark.verify(&air, &boundary, stark.deser_independent_proof_stream(proof.clone()));

            assert!(verdict.is_err(), "invalid stark proof verifies");
            println!("proof with invalid boundary rejected! \\o/");
        }

        // verify with false witness
        println!("attempting to prove with false witness (should fail) ...");

        let mut thread_rng = thread_rng();
        let mut rand_bytes = |n: usize| {
            let mut bytes = vec![0; n];
            thread_rng.fill_bytes(&mut bytes);
            bytes
        };

        let cycle = rand_bytes(1)[0] as usize % trace.len();
        let register = rand_bytes(1)[0] as usize % rp.m;
        let error = field.sample(&Bytes::new(rand_bytes(17)));

        trace[cycle][register] = trace[cycle][register] + error;

        let proof = stark.prove(trace, &air, &boundary, IndependentProofStream::new());
        assert!(proof.is_err(), "stark.prove should have failed with invalid trace");
    }
}
