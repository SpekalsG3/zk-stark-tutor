use crate::field::field::Field;
use crate::field::field_element::FieldElement;
use crate::field::polynomial::Polynomial;
use crate::m_polynomial::MPolynomial;
use crate::utils::matrix::{rref, transpose};

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
        m: usize, // 2
        _capacity: usize, // 1
        #[allow(non_snake_case)]
        N: usize, // 27
    ) -> Self {
        let g = field.smallest_generator();
        Self {
            field,
            m: 2,
            capacity: 1,
            N,
            alpha: g.value,
            alpha_inv: field.inv(field.neg_mod(g.value)), // alpha ^ (-1)
            MDS: Self::get_mds(g, m),
            MDS_inv: vec![
                vec![FieldElement::new(field, 210387253332845851216830350818816760948), FieldElement::new(field, 60110643809384528919094385948233360270)],
                vec![FieldElement::new(field, 90165965714076793378641578922350040407), FieldElement::new(field, 180331931428153586757283157844700080811)],
            ],
            // with_capacity = 2 * m * N
            round_constants: vec![174420698556543096520990950387834928928,109797589356993153279775383318666383471,228209559001143551442223248324541026000,268065703411175077628483247596226793933,250145786294793103303712876509736552288,154077925986488943960463842753819802236,204351119916823989032262966063401835731,57645879694647124999765652767459586992,102595110702094480597072290517349480965,8547439040206095323896524760274454544,50572190394727023982626065566525285390,87212354645973284136664042673979287772,64194686442324278631544434661927384193,23568247650578792137833165499572533289,264007385962234849237916966106429729444,227358300354534643391164539784212796168,179708233992972292788270914486717436725,102544935062767739638603684272741145148,65916940568893052493361867756647855734,144640159807528060664543800548526463356,58854991566939066418297427463486407598,144030533171309201969715569323510469388,264508722432906572066373216583268225708,22822825100935314666408731317941213728,33847779135505989201180138242500409760,146019284593100673590036640208621384175,51518045467620803302456472369449375741,73980612169525564135758195254813968438,31385101081646507577789564023348734881,270440021758749482599657914695597186347,185230877992845332344172234234093900282,210581925261995303483700331833844461519,233206235520000865382510460029939548462,178264060478215643105832556466392228683,69838834175855952450551936238929375468,75130152423898813192534713014890860884,59548275327570508231574439445023390415,43940979610564284967906719248029560342,95698099945510403318638730212513975543,77477281413246683919638580088082585351,206782304337497407273753387483545866988,141354674678885463410629926929791411677,19199940390616847185791261689448703536,177613618019817222931832611307175416361,267907751104005095811361156810067173120,33296937002574626161968730356414562829,63869971087730263431297345514089710163,200481282361858638356211874793723910968,69328322389827264175963301685224506573,239701591437699235962505536113880102063,17960711445525398132996203513667829940,219475635972825920849300179026969104558,230038611061931950901316413728344422823,149446814906994196814403811767389273580,25535582028106779796087284957910475912,93289417880348777872263904150910422367,4779480286211196984451238384230810357,208762241641328369347598009494500117007,34228805619823025763071411313049761059,158261639460060679368122984607245246072,65048656051037025727800046057154042857,134082885477766198947293095565706395050,23967684755547703714152865513907888630,8509910504689758897218307536423349149,232305018091414643115319608123377855094,170072389454430682177687789261779760420,62135161769871915508973643543011377095,15206455074148527786017895403501783555,201789266626211748844060539344508876901,179184798347291033565902633932801007181,9615415305648972863990712807943643216,95833504353120759807903032286346974132,181975981662825791627439958531194157276,267590267548392311337348990085222348350,49899900194200760923895805362651210299,89154519171560176870922732825690870368,265649728290587561988835145059696796797,140583850659111280842212115981043548773,266613908274746297875734026718148328473,236645120614796645424209995934912005038,265994065390091692951198742962775551587,59082836245981276360468435361137847418,26520064393601763202002257967586372271,108781692876845940775123575518154991932,138658034947980464912436420092172339656,45127926643030464660360100330441456786,210648707238405606524318597107528368459,42375307814689058540930810881506327698,237653383836912953043082350232373669114,236638771475482562810484106048928039069,168366677297979943348866069441526047857,195301262267610361172900534545341678525,2123819604855435621395010720102555908,96986567016099155020743003059932893278,248057324456138589201107100302767574618,198550227406618432920989444844179399959,177812676254201468976352471992022853250,211374136170376198628213577084029234846,105785712445518775732830634260671010540,122179368175793934687780753063673096166,126848216361173160497844444214866193172,22264167580742653700039698161547403113,234275908658634858929918842923795514466,189409811294589697028796856023159619258,75017033107075630953974011872571911999,144945344860351075586575129489570116296,261991152616933455169437121254310265934,18450316039330448878816627264054416127]
                .into_iter()
                .map(|v| FieldElement::new(field, v))
                .collect(),
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
        let rp = RescuePrime::new(&field, 2, 1, 27);

        assert_eq!(rp.alpha, 3);
        assert_eq!(rp.alpha_inv, 180331931428153586757283157844700080811);
    }

    #[test]
    fn hash() {
        let field = Field::new(FIELD_PRIME);
        let rp = RescuePrime::new(&field, 2, 1, 27);

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
        let rp = RescuePrime::new(&field, 2, 1, 27);

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
    fn constraints() {
        let field = Field::new(FIELD_PRIME);
        let rp = RescuePrime::new(&field, 2, 1, 27);

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
