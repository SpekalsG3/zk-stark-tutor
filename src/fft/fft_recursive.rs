use std::f32::consts::PI;
use crate::utils::complex::{Complex, I};

pub fn fft_recursive(
    inputs: Vec<Complex>,
) -> Vec<Complex> {
    fn fft_inner(inputs: Vec<Complex>) -> Vec<Complex> {
        let n = inputs.len();
        if n == 1 {
            // because:
            // we can't split n anymore and just compute y_k=\sum_{j=0}^{n-1}{...}
            // and n-1 == 0 means we have only one iteration with j=0 for x[j]*e^{j...}
            // which means x[0] * e^{0...} = x[0] * e^0 = x[0] * 1 = x[0]
            return inputs;
        }

        // divide
        let even = inputs
            .clone()
            .into_iter()
            .enumerate()
            .filter_map(|(i, x)| if i % 2 == 0 { Some(x) } else { None })
            .collect::<Vec<_>>();
        let odd  = inputs
            .into_iter()
            .enumerate()
            .filter_map(|(i, x)| if i % 2 == 1 { Some(x) } else { None })
            .collect::<Vec<_>>();

        let even = fft_inner(even);
        let odd  = fft_inner(odd);

        // conquer
        let mut bins = vec![Complex::zero(); n];
        for k in 0..n/2 {
            // having calculated the smallest coefficient (with n=1) we now have coefficients X_e and X_o
            // now we can use them `X_{k} = X_e[k] + o^k * X_o[k]`
            // and if we try to substitute k with k = k' + n/2 we will see that `^{+ n/2}` exponents cancel out
            // giving us X_{k'+h/2} = X_e[k'] - o^k * X_o[k']
            // then we can safely replace k' back to k: X_{k+h/2} = X_e[k] - o^k * X_o[k]
            // and we get the same expression but only with the minus instead

            let exp = I * -2_f32 * PI * (k as f32 / n as f32);
            let omega = exp.exp();
            let omega_x_odd = omega * odd[k].clone();
            bins[k] = even[k].clone() + omega_x_odd.clone();
            bins[k + n/2] = even[k].clone() - omega_x_odd;
        }

        bins
    }

    fft_inner(inputs)
}
