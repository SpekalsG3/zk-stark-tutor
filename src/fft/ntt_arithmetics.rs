use crate::fft::ntt::{intt, ntt};
use crate::field::field_element::FieldElement;
use crate::field::polynomial::Polynomial;

pub fn fast_multiply<'a>(
    root: FieldElement<'a>,
    root_order: usize,
    lhs: Polynomial<'a>,
    rhs: Polynomial<'a>,
) -> Polynomial<'a> {
    assert_eq!(
        root ^ root_order,
        root.field.one(),
        "supplied root {} does not have supplied root_order {}",
        root.value,
        root_order,
    );
    assert_ne!(
        root ^ (root_order / 2),
        root.field.one(),
        "supplied root {} is not a primitive of root_order {}",
        root.value,
        root_order,
    );

    if lhs.is_zero() || rhs.is_zero() {
        return Polynomial::new(vec![]);
    }

    let field = root.field;
    let mut root = root;
    let mut order = root_order;
    // SAFETY: checked if zero before
    let degree = lhs.degree().unwrap() + rhs.degree().unwrap();

    // for whatever reason `degree * 2 < order` doesnt work
    while degree < order / 2 {
        root = root ^ 2_usize;
        order = order / 2_usize;
    }

    let mut lhs = lhs.coefficients;
    lhs.extend((lhs.len()..order).map(|_| field.zero()));

    let mut rhs = rhs.coefficients;
    rhs.extend((rhs.len()..order).map(|_| field.zero()));

    let lhs = ntt(root, lhs);
    let rhs = ntt(root, rhs);

    let hadamard_product = (0..order)
        .map(|i| {
            // SAFETY: padded up to order elements
            lhs[i] * rhs[i]
        })
        .collect();

    let mut coeffs = intt(root, hadamard_product);
    if degree + 1 < coeffs.len() {
        coeffs.drain((degree + 1)..);
    }
    Polynomial::new(coeffs)
}

pub fn fast_zerofier<'a>(
    root: FieldElement<'a>,
    root_order: usize,
    domain: &[FieldElement<'a>],
) -> Polynomial<'a> {
    assert_eq!(
        root ^ root_order,
        root.field.one(),
        "supplied root {} does not have supplied root_order {}",
        root.value,
        root_order,
    );
    assert_ne!(
        root ^ (root_order / 2),
        root.field.one(),
        "supplied root {} is not a primitive of root_order {}",
        root.value,
        root_order,
    );

    fn inner<'i> (
        root: FieldElement<'i>,
        root_order: usize,
        domain: &[FieldElement<'i>],
    ) -> Polynomial<'i> {
        if domain.len() == 0 {
            return Polynomial::new(vec![]);
        }
        if domain.len() == 1 {
            return Polynomial::new(vec![
                -domain[0],
                root.field.one(),
            ]);
        }

        let half = domain.len() / 2;
        let left = inner(root, root_order, &domain[..half]);
        let right = inner(root, root_order, &domain[half..]);
        fast_multiply(root, root_order, left, right)
    }

    inner(root, root_order, domain)
}

pub fn fast_evaluate_domain<'a>(
    root: FieldElement<'a>,
    root_order: usize,
    polynomial: Polynomial<'a>,
    domain: &[FieldElement<'a>],
) -> Vec<FieldElement<'a>> {
    assert_eq!(
        root ^ root_order,
        root.field.one(),
        "supplied root {} does not have supplied root_order {}",
        root.value,
        root_order,
    );
    assert_ne!(
        root ^ (root_order / 2),
        root.field.one(),
        "supplied root {} is not a primitive of root_order {}",
        root.value,
        root_order,
    );

    fn inner<'i>(
        root: FieldElement<'i>,
        root_order: usize,
        polynomial: Polynomial<'i>,
        domain: &[FieldElement<'i>],
    ) -> Vec<FieldElement<'i>> {
        if domain.len() == 0 {
            return vec![];
        }
        if domain.len() == 1 {
            return vec![
                polynomial.evaluate(&domain[0]),
            ]
        }

        let half = domain.len() / 2;

        let left = fast_zerofier(root, root_order, &domain[..half]);
        let right = fast_zerofier(root, root_order, &domain[half..]);

        let mut left  = inner(root, root_order, polynomial.clone() % left, &domain[..half]);
        let     right = inner(root, root_order, polynomial % right, &domain[half..]);

        left.extend(right);
        left
    }

    inner(root, root_order, polynomial, domain)
}

pub fn fast_coset_evaluate<'a>(
    generator: FieldElement<'a>,
    root_order: usize,
    offset: FieldElement<'a>,
    polynomial: Polynomial<'a>,
) -> Vec<FieldElement<'a>> {
    let mut coeffs = polynomial.scale(offset).coefficients;
    coeffs.extend(vec![generator.field.zero(); root_order - coeffs.len()]);
    ntt(generator, coeffs)
}

pub fn fast_interpolate<'a>(
    root: FieldElement<'a>,
    root_order: usize,
    domain: &[FieldElement<'a>],
    values: &[FieldElement<'a>],
) -> Polynomial<'a> {
    assert_eq!(
        root ^ root_order,
        root.field.one(),
        "supplied root {} does not have supplied root_order {}",
        root.value,
        root_order,
    );
    assert_ne!(
        root ^ (root_order / 2),
        root.field.one(),
        "supplied root {} is not a primitive of root_order {}",
        root.value,
        root_order,
    );
    assert_eq!(domain.len(), values.len());

    fn inner<'i>(
        root: FieldElement<'i>,
        root_order: usize,
        domain: &[FieldElement<'i>],
        values: &[FieldElement<'i>],
    ) -> Polynomial<'i> {
        if domain.len() == 0 {
            return Polynomial::new(vec![]);
        }
        if domain.len() == 1 {
            return Polynomial::new(vec![values[0]]);
        }

        let half = domain.len() / 2;

        let left_zerofier = fast_zerofier(root, root_order, &domain[..half]);
        let right_zerofier = fast_zerofier(root, root_order, &domain[half..]);

        let left_offset  = fast_evaluate_domain(root, root_order, right_zerofier.clone(), &domain[..half]);
        let right_offset = fast_evaluate_domain(root, root_order, left_zerofier.clone(), &domain[half..]);

        let left_targets = left_offset
            .into_iter()
            .enumerate()
            .map(|(i, d)| {
                values[i] / d
            })
            .collect::<Vec<_>>();
        let right_targets = right_offset
            .into_iter()
            .enumerate()
            .map(|(i, d)| {
                values[i + half] / d
            })
            .collect::<Vec<_>>();

        let left_interpolant = inner(root, root_order, &domain[..half], &left_targets);
        let right_interpolant = inner(root, root_order, &domain[half..], &right_targets);

        left_interpolant * right_zerofier + right_interpolant * left_zerofier
    }

    inner(root, root_order, domain, values)
}

pub fn fast_coset_divide<'a>(
    root: FieldElement<'a>,
    root_order: usize,
    offset: FieldElement<'a>,
    lhs: Polynomial<'a>,
    rhs: Polynomial<'a>,
) -> Polynomial<'a> {
    assert_eq!(
        root ^ root_order,
        root.field.one(),
        "supplied root {} does not have supplied root_order {}",
        root.value,
        root_order,
    );
    assert_ne!(
        root ^ (root_order / 2),
        root.field.one(),
        "supplied root {} is not a primitive of root_order {}",
        root.value,
        root_order,
    );
    assert!(!rhs.is_zero(), "cannot divide by zero polynomial");

    if lhs.is_zero() {
        return Polynomial::new(vec![]);
    }

    // SAFETY: checked if zero before
    let lhs_degree = lhs.degree().unwrap();
    let rhs_degree = rhs.degree().unwrap();
    assert!(lhs_degree >= rhs_degree, "cannot divide by polynomial of larger degree");

    let field = root.field;
    let mut root = root;
    let mut order = root_order;
    let degree = lhs_degree.max(rhs_degree);
    let result_degree = lhs_degree - rhs_degree + 1;

    // for whatever reason `degree * 2 < order` doesnt work
    while degree < order / 2 {
        root = root ^ 2_usize;
        order = order / 2_usize;
    }

    let inner = |poly: Polynomial<'a>| {
        let poly = poly.scale(offset);
        let mut poly = poly.coefficients;
        poly.extend(
            (poly.len()..order)
                .map(|_| field.zero()),
        );
        ntt(root, poly)
    };

    let lhs = inner(lhs);
    let rhs = inner(rhs);

    let quotient = (0..order)
        .map(|i| {
            // SAFETY: padded up to order elements
            lhs[i] / rhs[i]
        })
        .collect();

    let mut coeffs = intt(root, quotient);
    if result_degree < coeffs.len() {
        coeffs.drain(result_degree..);
    }
    let scaled_poly = Polynomial::new(coeffs);

    scaled_poly.scale(offset.inverse())
}

#[cfg(test)]
mod tests {
    use rand::{RngCore, thread_rng};
    use rand::rngs::ThreadRng;
    use crate::fft::ntt_arithmetics::{fast_coset_divide, fast_coset_evaluate, fast_evaluate_domain, fast_interpolate, fast_multiply, fast_zerofier};
    use crate::field::field::{Field, FIELD_PRIME};
    use crate::field::field_element::FieldElement;
    use crate::field::polynomial::Polynomial;
    use crate::utils::bytes::Bytes;

    fn rand_domain<'a>(
        thread_rng: &mut ThreadRng,
        field: &'a Field,
        size: usize,
    ) -> Vec<FieldElement<'a>> {
        (0..size)
            .map(|_| {
                let mut bytes = [0; 17];
                thread_rng.fill_bytes(&mut bytes);
                field.sample(&Bytes::new(bytes.to_vec()))
            })
            .collect()
    }

    fn rand_poly<'a>(
        thread_rng: &mut ThreadRng,
        field: &'a Field,
        max_degree: usize,
    ) -> Polynomial<'a> {
        let degree = {
            let mut bytes = [0; 1];
            thread_rng.fill_bytes(&mut bytes);
            u8::from_be_bytes(bytes) as usize % max_degree
        };

        let coeffs = rand_domain(thread_rng, field, degree);
        Polynomial::new(coeffs)
    }

    #[test]
    fn multiply () {
        let field = Field::new(FIELD_PRIME);
        let n: usize = 1 << 6;
        let primitive_root = field.primitive_nth_root(n as u128);

        let mut thread_rng = thread_rng();

        for trial in 0..20 {
            let lhs = rand_poly(&mut thread_rng, &field, n / 2);
            let rhs = rand_poly(&mut thread_rng, &field, n / 2);

            assert_eq!(
                fast_multiply(primitive_root, n, lhs.clone(), rhs.clone()),
                lhs.clone() * rhs.clone(),
                "#{trial} trial failed with:\nlhs {:?}\nrhs {:?}",
                lhs.coefficients.iter().map(|f| f.value).collect::<Vec<_>>(),
                rhs.coefficients.iter().map(|f| f.value).collect::<Vec<_>>(),
            );
        }
    }

    #[test]
    fn zerofier () {
        let field = Field::new(FIELD_PRIME);
        let n: usize = 1 << 6;
        let primitive_root = field.primitive_nth_root(n as u128);

        let mut thread_rng = thread_rng();

        for trial in 0..20 {
            let poly = rand_poly(&mut thread_rng, &field, n);
            let zerofier = fast_zerofier(
                primitive_root,
                n,
                &poly.coefficients,
            );

            for c in &poly.coefficients {
                assert_eq!(
                    zerofier.evaluate(c),
                    field.zero(),
                    "#{trial} trial failed with:\npoly {:?}\nzerofier {:?}",
                    poly.coefficients.iter().map(|f| f.value).collect::<Vec<_>>(),
                    zerofier.coefficients.iter().map(|f| f.value).collect::<Vec<_>>(),
                )
            }
        }
    }

    #[test]
    fn evaluate_domain() {
        let field = Field::new(FIELD_PRIME);
        let n: usize = 1 << 6;
        let primitive_root = field.primitive_nth_root(n as u128);

        let mut thread_rng = thread_rng();

        for trial in 0..20 {
            let poly = rand_poly(&mut thread_rng, &field, n);
            let domain = rand_domain(&mut thread_rng, &field, n);

            let evaluation = fast_evaluate_domain(
                primitive_root,
                n,
                poly.clone(),
                &domain,
            );

            assert_eq!(
                poly.evaluate_domain(&domain),
                evaluation,
                "#{trial} trial failed with:\npoly {:?}\ndomain {:?}",
                poly.coefficients.iter().map(|f| f.value).collect::<Vec<_>>(),
                domain.iter().map(|f| f.value).collect::<Vec<_>>(),
            )
        }
    }

    #[test]
    fn interpolate() {
        let field = Field::new(FIELD_PRIME);
        let n: usize = 1 << 6;
        let primitive_root = field.primitive_nth_root(n as u128);

        let mut thread_rng = thread_rng();

        // let start = SystemTime::now();

        for trial in 0..20 {
            let domain = rand_domain(&mut thread_rng, &field, n);
            let values = rand_domain(&mut thread_rng, &field, n);

            // let poly = Polynomial::interpolate_domain(&domain, &values);
            // let evaluation = poly.evaluate_domain(&domain);
            let poly = fast_interpolate(primitive_root, n, &domain, &values);
            let evaluation = fast_evaluate_domain(
                primitive_root,
                n,
                poly.clone(),
                &domain,
            );

            assert_eq!(
                values,
                evaluation,
                "#{trial} trial failed with:\nvalues {:?}\ndomain {:?}",
                values.iter().map(|f| f.value).collect::<Vec<_>>(),
                domain.iter().map(|f| f.value).collect::<Vec<_>>(),
            )
        }

        // 31_809ms - default
        //  9_769ms - fast version
        // println!("done in {}ms", SystemTime::now().duration_since(start).unwrap().as_millis())
    }

    #[test]
    fn coset_evaluate () {
        let field = Field::new(FIELD_PRIME);
        let n: usize = 1 << 6;
        let primitive_root = field.primitive_nth_root(n as u128);

        let offset = FieldElement::new(&field, 5);
        let domain = (0..n)
            .map(|i| (primitive_root ^ i) * offset)
            .collect::<Vec<_>>();

        let mut thread_rng = thread_rng();
        let poly = rand_poly(&mut thread_rng, &field, n);

        // 115ms
        let evaluation_domain = fast_evaluate_domain(primitive_root, n, poly.clone(), &domain);
        // 1ms
        let evaluation_coset = fast_coset_evaluate(primitive_root, n, offset, poly);

        assert_eq!(evaluation_domain, evaluation_coset);
    }

    #[test]
    fn coset_divide() {
        let field = Field::new(FIELD_PRIME);
        let n: usize = 1 << 6;
        let primitive_root = field.primitive_nth_root(n as u128);

        let mut thread_rng = thread_rng();

        for trial in 0..20 {
            let lhs = rand_poly(&mut thread_rng, &field, n / 2);
            let rhs = rand_poly(&mut thread_rng, &field, n / 2);

            let prod = fast_multiply(primitive_root, n, lhs.clone(), rhs.clone());
            let div = fast_coset_divide(primitive_root, n, field.generator(), prod, lhs.clone());

            assert_eq!(
                div,
                rhs,
                "#{trial} trial failed with:\nlhs {:?}\nrhs {:?}",
                lhs.coefficients.iter().map(|f| f.value).collect::<Vec<_>>(),
                rhs.coefficients.iter().map(|f| f.value).collect::<Vec<_>>(),
            )
        }
    }
}
