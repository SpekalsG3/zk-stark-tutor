pub fn gcd (a: u128, b: u128) -> u128 {
  let mut a = a;
  let mut b = b;
  while b != 0 {
    let t = b;
    b = a % b;
    a = t;
  }
  a
}

#[cfg(test)]
mod tests {
  use crate::field::field::FIELD_PRIME;
  use crate::utils::gcd::gcd;

  #[test]
  fn test () {
    assert_eq!(gcd(10, 5), 5);
    assert_eq!(gcd(240, 46), 2);
    assert_eq!(gcd(3, FIELD_PRIME), 1);
  }
}
