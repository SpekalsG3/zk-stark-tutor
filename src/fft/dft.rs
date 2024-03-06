use std::f32::consts::PI;
use crate::utils::complex::{Complex, I};

pub fn dft (inputs: &[f32]) -> Vec<Complex> {
    let n = inputs.len();

    (0..n )
        .map(|f| {
            inputs
                .iter()
                .enumerate()
                .map(|(i, x)| {
                    // minus is to calculate fourier coefficient
                    // without minus it's an inverse operation - inverse DFT
                    // but i have no idea why result is the same xd
                    let exp = I * -2_f32 * PI * i as f32 * (f as f32 / n as f32);
                    let omega = exp.exp();
                    omega * *x
                })
                .reduce(|a, b| a + b)
                .unwrap()
        })
        .collect()
}
