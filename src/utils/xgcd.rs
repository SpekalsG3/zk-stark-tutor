use std::ops::Neg;

// original from tutorial
#[deprecated(note="overflow when working close to boundary, use `u_xgcd`")]
pub fn xgcd (a: u128, b: u128) -> (i128, i128, u128) {
  let mut r = (a, b);
  let mut s = (1_i128, 0_i128);
  let mut t = (0_i128, 1_i128);

  while r.1 != 0 {
    let q = r.0 / r.1;
    r = (r.1, r.0 - q * r.1);
    s = (s.1, s.0 - q as i128 * s.1);
    t = (t.1, t.0 - q as i128 * t.1);
  }

  (s.0, t.0, r.0) // a, b, g
}

// version using unsigned arithmetics - no overflow
// source - https://jeffhurchalla.com/2018/10/13/implementing-the-extended-euclidean-algorithm-with-unsigned-inputs/
pub fn u_xgcd (a: u128, b: u128) -> (i128, i128, u128) {
  let mut xy1 = (1_i128, 0_i128);

  let mut xy0 = (0_i128, 1_i128);

  let mut a = (a, b);

  let mut q = 0_u128;

  while a.1 != 0 {
    let x2 = xy0.0 - (q as i128) * xy1.0;
    let y2 = xy0.1 - (q as i128) * xy1.1;

    xy0 = xy1;
    xy1 = (x2, y2);

    // let a0 = a.0;
    // a.0 = a.1;
    // q = a0 / a.0;
    // a.1 = a0 - q * a.0;

    q = a.0 / a.1;
    a = (a.1, a.0 - q * a.1)
  }

  (xy1.0, xy1.1, a.0)
}

pub fn multiplicative_inverse (x: u128, y: u128, m: u128) -> u128 {
  let (a, _, _) = u_xgcd(y, m);
  if a < 0 {
    (m - x) * (a.neg() as u128) % m
  } else {
    x * (a as u128) % m
  }
}

#[cfg(test)]
mod tests {
  use crate::field::field::FIELD_PRIME;
  use super::*;

  #[test]
  #[allow(deprecated)]
  fn test_xgcd () {
    assert_eq!(xgcd(10, 5), (0, 1, 5));
    assert_eq!(xgcd(240, 46), (-9, 47, 2));
    // assert_eq!(xgcd(3, 1 + BIG_PRIME), (1, 0, 1)); // fails
  }

  #[test]
  fn test_u_xgcd () {
    assert_eq!(u_xgcd(10, 5), (0, 1, 5));
    assert_eq!(u_xgcd(240, 46), (-9, 47, 2));
    assert_eq!(u_xgcd(3, FIELD_PRIME), (90165965714076793378641578922350040406, -1, 1));
  }

  #[test]
  fn u_xgcd_multiplicative_inverse () {
    assert_eq!(multiplicative_inverse( 6, 2, 23), 3);
    assert_eq!(multiplicative_inverse(12, 4, 23), 3);
    assert_eq!(multiplicative_inverse( 8, 8, 23), 1);
    assert_eq!(multiplicative_inverse( 1, 1, 23), 1);
    assert_eq!(multiplicative_inverse( 0, 5, 23), 0);
  }
}
