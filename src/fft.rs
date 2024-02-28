use std::f32::consts::PI;
use std::ops::{Add, Mul, Sub};
use crate::utils::bit_iter::BitIter;

#[derive(Debug, Clone, PartialEq)]
pub struct Complex {
    re: f32,
    im: f32,
}

impl Complex {
    pub fn zero () -> Self {
        Self::new(0_f32, 0_f32)
    }
    pub fn one () -> Self {
        Self::new(1_f32, 0_f32)
    }
    pub fn new (re: f32, im: f32) -> Self {
        Self { re, im }
    }

    pub fn from_polar(r: f32, theta: f32) -> Self {
        Self::new(r * theta.cos(), r * theta.sin())
    }

    /// Computes `e^(self)`, where `e` is the base of the natural logarithm.
    #[inline]
    pub fn exp(self) -> Self {
        // formula: e^(a + bi) = e^a (cos(b) + i*sin(b)) = from_polar(e^a, b)

        let Complex { re, mut im } = self;
        // Treat the corner cases +∞, -∞, and NaN
        if re.is_infinite() {
            if re < 0_f32 {
                if !im.is_finite() {
                    return Self::new(0_f32, 0_f32);
                }
            } else {
                if im == 0_f32 || !im.is_finite() {
                    if im.is_infinite() {
                        im = f32::NAN;
                    }
                    return Self::new(re, im);
                }
            }
        } else if re.is_nan() && im == 0_f32 {
            return self;
        }

        Self::from_polar(re.exp(), im)
    }
}
impl Mul<Complex> for Complex {
    type Output = Self;

    #[inline]
    fn mul(self, other: Self) -> Self::Output {
        let re = self.re * other.re - self.im * other.im;
        let im = self.re * other.im + self.im * other.re;
        Self::Output::new(re, im)
    }
}
// (a + i b) + (c + i d) == (a + c) + i (b + d)
impl Add<Complex> for Complex {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self::Output::new(self.re + other.re, self.im + other.im)
    }
}
// (a + i b) + (c + i d) == (a + c) + i (b + d)
impl Sub<Complex> for Complex {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self::Output::new(self.re - other.re, self.im - other.im)
    }
}
impl Mul<f32> for Complex {
    type Output = Complex;

    #[inline]
    fn mul(self, other: f32) -> Self::Output {
        Self::Output::new(self.re * other, self.im * other)
    }
}

const I: Complex = Complex { re: 0_f32, im: 1_f32 };

pub fn fft_recursive(inputs: &[f32]) -> Vec<Complex> {
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

    fft_inner(inputs.into_iter().map(|x| Complex::new(*x, 0_f32)).collect())
}

pub fn fft(inputs: Vec<Complex>) -> Vec<Complex> {
    fn fft_inner(inputs: &mut Vec<Complex>) {
        let n = inputs.len();

        let mut m = 1;
        while m < n {
            m = m << 1;

            let o_k = {
                let exp = I * -2_f32 * (PI / m as f32);
                exp.exp()
            };
            for k in (0..n).step_by(m) {
                let mut o = Complex::one();

                for j in 0..m/2 {
                    let odd = o.clone() * inputs[k + j + m/2].clone();
                    let even = inputs[k + j].clone();
                    inputs[k + j] = even.clone() + odd.clone();
                    inputs[k + j + m/2] = even - odd;
                    o = o * o_k.clone();
                }
            }
        }
    }

    let orig_n = inputs.len();
    let n = orig_n.next_power_of_two();
    let mut new_inputs = vec![Complex::zero(); n];

    // pad the input with zeros if the input is not the size of the power two
    let inputs = inputs
        .into_iter()
        .chain(
            vec![Complex::zero(); n - orig_n]
        );

    // bit-reverse-copy
    let biggest_bit = (n - 1).ilog2() as usize;
    for (k, el) in inputs.enumerate() {
        let k_rev = {
            let mut n = 0;
            let bit_iter = BitIter(Some(biggest_bit), k);
            for (i, bit) in bit_iter.enumerate() {
                n = n + (if bit { 1 } else { 0 } << i)
            }
            n
        };
        new_inputs[k_rev] = el;
    }

    fft_inner(&mut new_inputs);
    new_inputs
}

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

#[cfg(test)]
mod tests {
    // use std::f32::consts::PI;
    use std::time::SystemTime;
    use crate::fft::{
        Complex,
        // dft,
        fft,
        // fft_recursive,
    };

    #[test]
    fn test () {
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
        // let freqs = fft_recursive(&inputs);
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

