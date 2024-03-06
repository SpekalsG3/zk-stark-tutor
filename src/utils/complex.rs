use std::ops::{Add, Mul, Sub};

#[derive(Debug, Clone, PartialEq)]
pub struct Complex {
    pub re: f32,
    pub im: f32,
}

pub const I: Complex = Complex { re: 0_f32, im: 1_f32 };

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
