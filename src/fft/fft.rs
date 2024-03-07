use std::f32::consts::PI;
use crate::utils::bit_reverse_copy::bit_reverse_copy;
use crate::utils::complex::{Complex, I};

pub fn fft(inputs: Vec<Complex>) -> Vec<Complex> {
    let mut inputs = bit_reverse_copy(Complex::zero(), inputs);
    let n = inputs.len();

    // iterative fft
    let mut size = 1;
    while size < n {
        size = size << 1;

        let halfsize = size / 2;
        let o_k = {
            let exp = I * -2_f32 * (PI / size as f32);
            exp.exp()
        };
        for k in (0..n).step_by(size) {
            let mut o = Complex::one();

            for j in k..k+halfsize {
                let odd = o.clone() * inputs[j + halfsize].clone();
                let even = inputs[j].clone();
                inputs[j           ] = even.clone() + odd.clone();
                inputs[j + halfsize] = even - odd;
                o = o * o_k.clone();
            }
        }
    }

    inputs
}

#[cfg(test)]
mod tests {
    use std::time::SystemTime;
    use crate::fft::fft::fft;
    use crate::utils::complex::Complex;

    #[test]
    fn test_ffts () {
        let n: usize = 8;
        println!("{n} samples");
        let inputs = vec![1_f32, 1_f32, 1_f32, 1_f32, 0_f32, 0_f32, 0_f32, 0_f32]
            .into_iter()
            .map(|x| Complex::new(x, 0_f32))
            .collect();
        // let inputs = (0..n)
        //     .map(|i| {
        //         let c = 2.0 * PI * (i as f32 / n as f32);
        //         (c * 1.0).cos() + (c * 2.0).sin()
        //     })
        //     .collect::<Vec<_>>();

        // inputs.iter().enumerate().for_each(|(i, x)| { println!("{:2}: {:.3}", i, x) });

        // let start = SystemTime::now();
        // let freqs = dft(&inputs);
        // println!("dft ({}ms):", SystemTime::now().duration_since(start).unwrap().as_millis()); // 1121ms for 4096 samples
        // freqs.iter().enumerate().for_each(|(i, x)| { println!("{:2}: ({:.3}, {:.3})", i, x.re, x.im) });

        // let start = SystemTime::now();
        // let freqs = fft_recursive(inputs);
        // println!("fft_recursive ({}ms):", SystemTime::now().duration_since(start).unwrap().as_millis()); // 8ms for 4096 samples
        // freqs.iter().enumerate().for_each(|(i, x)| { println!("{:2}: ({:.3}, {:.3})", i, x.re, x.im) });

        let start = SystemTime::now();
        let freqs = fft(inputs);
        println!("fft ({}ms):", SystemTime::now().duration_since(start).unwrap().as_millis()); // 1ms for 4096 samples
        // freqs.iter().enumerate().for_each(|(i, x)| { println!("{:2}: ({:.3}, {:.3})", i, x.re, x.im) });

        let correct = vec![
            Complex::new(4.000,  0.000),
            Complex::new(1.000, -2.414),
            Complex::new(0.000,  0.000),
            Complex::new(1.000, -0.414),
            Complex::new(0.000,  0.000),
            Complex::new(1.000,  0.414),
            Complex::new(0.000,  0.000),
            Complex::new(1.000,  2.414),
        ];
        for i in 0..n {
            // compare strings because float values are not exact
            let a = format!("({:.3}, {:.3})", freqs[i].re, freqs[i].im);
            let b = format!("({:.3}, {:.3})", correct[i].re, correct[i].im);
            assert_eq!(a, b);
        }
    }
}
