use std::ops::Add;

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct Bytes(Vec<u8>);

impl Bytes {
  pub fn new (buf: Vec<u8>) -> Self {
    Self(buf)
  }

  pub fn to_hex (&self) -> String {
    self.0
      .iter()
      .map(|b| format!("{:02x}", b))
      .collect::<String>()
  }

  pub fn bytes (&self) -> &[u8] {
    &self.0
  }

  pub fn stringify (&self) -> String {
    format!("\"{}\"", self.to_hex())
  }

  pub fn unstringify (str: &str) -> Self {
    Self::from(&str[1..str.len()-1])
  }
}

impl From<&str> for Bytes {
  fn from (s: &str) -> Self {
    Bytes (
      (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
        .collect()
    )
  }
}

impl From<&[u8]> for Bytes {
  fn from (value: &[u8]) -> Self {
    Bytes (
      value.to_vec(),
    )
  }
}

impl From<Vec<u8>> for Bytes {
  fn from(value: Vec<u8>) -> Self {
    Self(value)
  }
}

impl Add for Bytes {
  type Output = Bytes;
  fn add (mut self, rhs: Self) -> Self::Output {
    self.0.extend(rhs.0);
    self
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn to_hex () {
    let bytes = Bytes(vec![0x49,0x6e,0x20,0x74]);
    assert_eq!(bytes.to_hex(), String::from("496e2074"));
  }

  #[test]
  fn from_str () {
    let bytes = Bytes (vec![0x49,0x6e,0x20,0x74]);
    assert_eq!(bytes, "496e2074".into());
  }

  #[test]
  fn add () {
    let bytes_a = Bytes (vec![0x49,0x6e]);
    let bytes_b = Bytes (vec![0x20,0x74]);
    assert_eq!(bytes_a + bytes_b, Bytes (vec![0x49,0x6e,0x20,0x74]));
  }
}
