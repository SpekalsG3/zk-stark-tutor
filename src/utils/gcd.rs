pub fn gcd (a: usize, b: usize) -> usize {
  let mut a = a;
  let mut b = b;
  while b != 0 {
    let t = b;
    b = a % b;
    a = t;
  }
  a
}
