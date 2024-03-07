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

#[cfg(test)]
mod tests {
    use rand::{RngCore, thread_rng};
    use rand::rngs::ThreadRng;
    use crate::fft::ntt_arithmetics::{fast_multiply};
    use crate::field::field::{Field, FIELD_PRIME};
    use crate::field::polynomial::Polynomial;
    use crate::utils::bytes::Bytes;

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

        let coefs = (0..degree)
            .map(|_| {
                let mut bytes = [0; 17];
                thread_rng.fill_bytes(&mut bytes);
                field.sample(&Bytes::new(bytes.to_vec()))
            })
            .collect();
        Polynomial::new(coefs)
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
}
