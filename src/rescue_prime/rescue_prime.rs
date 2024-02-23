use crate::crypto::shake256::shake256;
use crate::field::field::Field;
use crate::field::field_element::FieldElement;
use crate::field::polynomial::Polynomial;
use crate::m_polynomial::MPolynomial;
use crate::utils::bit_iter::BitIter;
use crate::utils::bytes::Bytes;
use crate::utils::matrix::{inverse, rref, transpose};

#[allow(non_snake_case)]
pub struct RescuePrime<'a> {
    field: &'a Field,
    pub(crate) m: usize,
    capacity: usize,
    pub(crate) N: usize,
    alpha: u128,
    alpha_inv: u128,
    MDS: Vec<Vec<FieldElement<'a>>>,
    MDS_inv: Vec<Vec<FieldElement<'a>>>,
    round_constants: Vec<FieldElement<'a>>,
}

struct IterPermutation<'a> {
    rescue_prime: &'a RescuePrime<'a>,
    pub state: Vec<FieldElement<'a>>,
    r: usize,
}

impl<'a> IterPermutation<'a> {
    pub fn new(
        rescue_prime: &'a RescuePrime<'a>,
        state: Vec<FieldElement<'a>>,
    ) -> Self {
        assert_eq!(state.len(), rescue_prime.capacity, "Received wrong number of input elements");

        // absorb
        let mut state = state;
        state.extend(
            vec![rescue_prime.field.zero(); rescue_prime.m - rescue_prime.capacity]
        );

        Self {
            rescue_prime,
            state,
            r: 0,
        }
    }
}

impl<'a> Iterator for IterPermutation<'a> {
    type Item = ();
    fn next(&mut self) -> Option<Self::Item> {
        if self.r == self.rescue_prime.N {
            return None;
        }

        self.state = self.state
            .iter()
            .enumerate()
            // S-box - alpha exponent
            .map(|(i, s)| {
                (i, s.clone() ^ self.rescue_prime.alpha)
            })
            // matrix
            .fold(
                vec![self.rescue_prime.field.zero(); self.rescue_prime.m],
                |mut acc, (i, s)| {
                    for j in 0..self.rescue_prime.m {
                        acc[j] = acc[j] + self.rescue_prime.MDS[j][i] * s
                    }
                    acc
                },
            )
            .into_iter()
            .enumerate()
            // constants
            .map(|(i, s)| {
                (i, s + self.rescue_prime.round_constants[2 * self.r * self.rescue_prime.m + i])
            })
            // inverse S-box - inverse alpha exponent
            .map(|(i, s)| (i, s ^ self.rescue_prime.alpha_inv))
            // matrix
            .fold(
                vec![self.rescue_prime.field.zero(); self.rescue_prime.m],
                |mut acc, (i, s)| {
                    for j in 0..self.rescue_prime.m {
                        acc[j] = acc[j] + self.rescue_prime.MDS[j][i] * s
                    }
                    acc
                },
            )
            .into_iter()
            .enumerate()
            // constants
            .map(|(i, s)| {
                s + self.rescue_prime.round_constants[2 * self.r * self.rescue_prime.m + self.rescue_prime.m + i]
            })
            .collect::<Vec<_>>();

        self.r += 1;

        Some(())
    }
}

impl<'a> RescuePrime<'a> {
    pub fn new (
        field: &'a Field,
        m: usize,
        capacity: usize,
        #[allow(non_snake_case)]
        N: usize,
        security_level: usize,
    ) -> Self {
        let g = field.smallest_generator();
        let mds = Self::get_mds(g, m);
        Self {
            field,
            m: 2,
            capacity,
            N,
            alpha: g.value,
            alpha_inv: field.inv(field.neg_mod(g.value)), // alpha ^ (-1)
            MDS: mds.clone(),
            MDS_inv: inverse(field, mds).expect("MDS should be invertible"),
            round_constants: Self::get_round_constants(field, m, capacity, security_level, N),
        }
    }

    fn get_mds<'m>(g: FieldElement<'m>, m: usize) -> Vec<Vec<FieldElement<'m>>> {
        let mut matrix = (0..m)
            .map(|i| {
                (0..(2 * m))
                    .map(|j| g ^ (i * j))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        rref(&mut matrix);

        let matrix = matrix
            .into_iter()
            .map(|row| {
                row[m..].to_vec()
            })
            .collect::<Vec<_>>();

        let matrix = transpose(matrix);

        return matrix
    }

    fn get_round_constants(field: &Field, m: usize, capacity: usize, security_level: usize, N: usize) -> Vec<FieldElement<'_>> {
        // generate pseudorandom bytes
        let bytes_per_int = (BitIter::from(field.order).count() + 7) / 8 + 1; // `+7` => ceil
        let num_bytes = bytes_per_int * 2 * m * N;
        let seed_string = format!("Rescue-XLIX({},{},{},{})", field.order, m, capacity, security_level);
        let bytes = shake256(Bytes::from(seed_string.as_bytes()), num_bytes);
        let bytes = bytes.bytes();

        // process byte string in chunks
        let f256 = FieldElement::new(field, 256);
        (0..2*m*N)
            .map(|i| {
                let chunk = &bytes[ (bytes_per_int * i)..(bytes_per_int * (i+1)) ];
                chunk
                    .into_iter()
                    .enumerate()
                    .map(|(j, b)| {
                        (f256 ^ j) * FieldElement::new(field, *b as u128)
                    })
                    .reduce(|a, b| a + b)
                    .unwrap()
            })
            .collect()
    }


    pub fn hash<'m> (&'m self, input_elements: FieldElement<'m>) -> FieldElement<'m> {
        let mut iter = IterPermutation::new(self, vec![input_elements]);
        iter.by_ref().last().unwrap();

        // squeeze
        iter.state.truncate(self.capacity);
        iter.state[0]
    }

    pub fn trace<'m> (&'m self, input_elements: FieldElement<'m>) -> Vec<Vec<FieldElement<'m>>> {
        let mut iter = IterPermutation::new(self, vec![input_elements]);

        let mut vec = Vec::with_capacity(self.N + 1);
        vec.push(iter.state.clone());

        while let Some(_) = iter.next() {
            vec.push(iter.state.clone())
        }
        vec
    }

    fn round_constants_polynomials<'m> (&'m self, omicron: FieldElement<'m>)
        -> (Vec<MPolynomial<'m>>, Vec<MPolynomial<'m>>) {
        let domain = (0..self.N)
            .map(|r| omicron ^ r)
            .collect::<Vec<_>>();

        let left = (0..self.m)
            .map(|i| {
                let values = (0..self.N)
                    .map(|r| self.round_constants[2 * r * self.m + i])
                    .collect::<Vec<_>>();
                let poly = Polynomial::interpolate_domain(&domain, &values);
                MPolynomial::lift(&poly, 0)
            })
            .collect();

        let right = (0..self.m)
            .map(|i| {
                let values = (0..self.N)
                    .map(|r| self.round_constants[2 * r * self.m + self.m + i])
                    .collect::<Vec<_>>();
                let poly = Polynomial::interpolate_domain(&domain, &values);
                MPolynomial::lift(&poly, 0)
            })
            .collect();

        (left, right)
    }

    pub fn transition_constraints<'m> (&'m self, omicron: FieldElement<'a>) -> Vec<MPolynomial<'m>> {
        // get polynomials that interpolate through the round constants
        let (first_step, second_step) = self.round_constants_polynomials(omicron);

        // arithmetize one round of Rescue-Prime
        let variables = MPolynomial::variables(1 + 2 * self.m, self.field);
        let previous_state = &variables[1..(1 + self.m)];
        let next_state = &variables[(1 + self.m)..(1 + 2 * self.m)];

        (0..self.m)
            .map(|i| {
                // compute left hand side symbolically
                let lhs = (0..self.m)
                    .map(|k| {
                        MPolynomial::constant(self.MDS[i][k]) * (previous_state[k].clone() ^ self.alpha)
                    })
                    .reduce(|a, b| a + b)
                    .unwrap() + first_step[i].clone();

                let rhs = (0..self.m)
                    .map(|k| {
                        MPolynomial::constant(self.MDS_inv[i][k]) * (next_state[k].clone() - second_step[k].clone())
                    })
                    .reduce(|a, b| a + b)
                    .unwrap() ^ self.alpha;

                lhs - rhs
            })
            .collect()
    }

    pub fn boundary_constraints<'m> (&'m self, output_element: FieldElement<'m>) -> Vec<(usize, usize, FieldElement<'m>)> {
        vec![
            (0, 1, self.field.zero()), // at start, capacity is zero
            (self.N, 0, output_element), // at end, rate part is the given output element
        ]
    }
}

#[cfg(test)]
mod tests {
    use rand::{RngCore, thread_rng};
    use crate::field::field::{Field, FIELD_PRIME};
    use crate::field::field_element::FieldElement;
    use crate::rescue_prime::rescue_prime::RescuePrime;
    use crate::utils::bytes::Bytes;

    #[test]
    fn new() {
        let field = Field::new(FIELD_PRIME);
        let rp = RescuePrime::new(&field, 2, 1, 128, 27);

        assert_eq!(rp.alpha, 3);
        assert_eq!(rp.alpha_inv, 180331931428153586757283157844700080811);
    }

    #[test]
    fn hash() {
        let field = Field::new(FIELD_PRIME);
        let rp = RescuePrime::new(&field, 2, 1, 128, 27);

        assert_eq!(
            rp.hash(FieldElement::new(&field, 1)),
            FieldElement::new(&field, 244180265933090377212304188905974087294)
        );
        assert_eq!(
            rp.hash(FieldElement::new(&field, 1)),
            FieldElement::new(&field, 244180265933090377212304188905974087294)
        );
    }

    #[test]
    fn trace() {
        let field = Field::new(FIELD_PRIME);
        let rp = RescuePrime::new(&field, 2, 1, 128, 27);

        let a = FieldElement::new(&field, 57322816861100832358702415967512842988);
        let b = FieldElement::new(&field, 89633745865384635541695204788332415101);
        let trace = rp.trace(a);
        assert!(trace[0][0] == a && trace[trace.len() - 1][0] == b, "rescue prime trace does not satisfy boundary conditions");
    }

    #[test]
    fn get_mds () {
        let field = Field::new(FIELD_PRIME);
        let mds = RescuePrime::get_mds(field.smallest_generator(), 2);
        assert_eq!(mds, vec![
            vec![FieldElement::new(&field, 270497897142230380135924736767050121214), FieldElement::new(&field,  4)],
            vec![FieldElement::new(&field, 270497897142230380135924736767050121205), FieldElement::new(&field, 13)],
        ])
    }

    #[test]
    fn get_round_constants() {
        let field = Field::new(FIELD_PRIME);
        let constants = RescuePrime::get_round_constants(&field, 2, 1, 128, 27);
        println!("{:?}", constants);
    }

    #[test]
    fn constraints() {
        let field = Field::new(FIELD_PRIME);
        let rp = RescuePrime::new(&field, 2, 1, 128, 27);

        let input = FieldElement::new(&field, 57322816861100832358702415967512842988);
        let output = rp.hash(input.clone());
        assert_eq!(output, FieldElement::new(&field, 89633745865384635541695204788332415101));

        let mut trace = rp.trace(input.clone());

        #[derive(Debug, PartialEq)]
        enum CheckResult {
            BoundaryError,
            TransitionError,
            Ok,
        }
        let check_constraints = |trace: &Vec<Vec<FieldElement>>| -> CheckResult {
            for (cycle, element, value) in rp.boundary_constraints(output) {
                if trace[cycle][element] != value {
                    return CheckResult::BoundaryError;
                }
            }

            let omicron = field.primitive_nth_root(1 << 119);
            for i in 0..(trace.len() - &1) {
                for poly in rp.transition_constraints(omicron) {
                    let mut point = Vec::with_capacity(1 + 2 * 2 * rp.m);
                    point.push(omicron ^ i);
                    point.extend(&trace[i]); // prev state
                    point.extend(&trace[i + 1]); // next state

                    if poly.evaluate(&point) != field.zero() {
                        return CheckResult::TransitionError;
                    }
                }
            }

            CheckResult::Ok
        };

        assert_eq!(check_constraints(&trace), CheckResult::Ok);

        println!("test invalid trace");

        let mut thread_rng = thread_rng();
        let mut rand_bytes = |n: usize| {
            let mut bytes = vec![0; n];
            thread_rng.fill_bytes(&mut bytes);
            bytes
        };

        let mut tests = Vec::with_capacity(10);
        tests.push((22, 1, FieldElement::new(&field, 17274817952119230544216945715808633996)));
        tests.resize_with(10, || {
            loop {
                let cycle = rand_bytes(1)[0] as usize % (rp.N + 1);
                let register = rand_bytes(1)[0] as usize % rp.m;
                let value = field.sample(&Bytes::new(rand_bytes(17)));

                if value.is_zero() { // because zero doesn't affect the trace
                    continue
                }

                return (cycle, register, value)
            }
        });

        for (i, (cycle, register, value)) in tests.into_iter().enumerate() {
            println!("test #{i}");

            trace[cycle][register] = trace[cycle][register] + value;

            assert_ne!(check_constraints(&trace), CheckResult::Ok, "\ncycle:\t{}\nregister:\t{}\nvalue:\t{}\n", cycle, register, value.value);

            trace[cycle][register] = trace[cycle][register] - value; // reset back to valid trace
        }
    }
}
