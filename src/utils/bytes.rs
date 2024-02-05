use std::ops::Add;
use serde::{Serialize, Serializer};

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct Bytes {
  bytes: Vec<u8>
}

impl Serialize for Bytes {
  // todo tests
  fn serialize<S> (&self, serializer: S) -> Result<S::Ok, S::Error> where S: Serializer {
    self.to_hex().serialize(serializer)
  }
}

pub struct Iter<'a> {
  i: usize,
  bytes: &'a Bytes,
}

impl Bytes {
  pub fn new (buf: Vec<u8>) -> Self {
    Self {
      bytes: buf.to_vec()
    }
  }

  pub fn to_hex (&self) -> String {
    self
      .bytes
      .iter()
      .map(|b| format!("{:02x}", b))
      .collect::<String>()
  }

  pub fn bytes (&self) -> &[u8] {
    &self.bytes
  }

  pub fn iter (&self) -> Iter<'_> {
    Iter {
      i: 0,
      bytes: self,
    }
  }
}

impl From<&str> for Bytes {
  fn from (s: &str) -> Self {
    Bytes {
      bytes: (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
        .collect(),
    }
  }
}

impl From<&[u8]> for Bytes {
  fn from (value: &[u8]) -> Self {
    Bytes {
      bytes: value.to_vec(),
    }
  }
}

impl Add for Bytes {
  type Output = Bytes;
  fn add (mut self, rhs: Self) -> Self::Output {
    self.bytes.extend(rhs.bytes);
    self
  }
}

impl<'a> Iterator for Iter<'a> {
  type Item = &'a u8;
  fn next (&mut self) -> Option<Self::Item> {
    let item = self.bytes.bytes.get(self.i);

    self.i += 1;

    item
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_to_hex () {
    let bytes = Bytes {
      bytes: vec![0x49,0x6e,0x20,0x74],
    };
    assert_eq!(bytes.to_hex(), String::from("496e2074"));
  }

  #[test]
  fn test_from_str () {
    let bytes = Bytes {
      bytes: vec![0x49,0x6e,0x20,0x74],
    };
    assert_eq!(bytes, "496e2074".into());
  }

  #[test]
  fn add () {
    let bytes_a = Bytes {
      bytes: vec![0x49,0x6e],
    };
    let bytes_b = Bytes {
      bytes: vec![0x20,0x74],
    };
    assert_eq!(bytes_a + bytes_b, Bytes {
      bytes: vec![0x49,0x6e,0x20,0x74],
    });
  }

  #[test]
  fn iter () {
    let bytes = Bytes {
      bytes: vec![0x49,0x6e,0x20,0x74],
    };
    let mut iter = bytes.iter();
    assert_eq!(iter.next(), Some(&0x49));
    assert_eq!(iter.next(), Some(&0x6e));
    assert_eq!(iter.next(), Some(&0x20));
    assert_eq!(iter.next(), Some(&0x74));
    assert_eq!(iter.next(), None);
  }
}
