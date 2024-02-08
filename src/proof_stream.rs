use std::fmt::{Display, Formatter};
use std::io::Read;
use serde::{Serialize};
use sha3::{digest::{Update, ExtendableOutput}, Shake256};
use crate::utils::bytes::Bytes;

#[derive(Debug)]
pub struct ProofStream<T> {
  objects: Vec<T>,
  read_index: usize,
}

impl<T> Display for ProofStream<T>
  where
    T: Display
{
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    let mut iter = self.objects.iter();
    if let Some(o) = iter.next() {
      write!(f, "{}", o)?;
    }
    for o in iter {
      write!(f, ",{}", o)?;
    }
    Ok(())
  }
}

impl<T> ProofStream<T>
  where
    T: Serialize + Clone
{
  pub fn new () -> Self {
    ProofStream {
      objects: vec![],
      read_index: 0,
    }
  }

  pub fn push (&mut self, obj: T) {
    self.objects.push(obj);
  }

  pub fn pull (&mut self) -> Option<T> {
    assert!(self.read_index < self.objects.len(), "Cannot pull, queue is empty");

    let obj = match self.objects.get(self.read_index) {
      Some(o) => Some(o.clone()),
      None => None
    };
    self.read_index += 1;

    obj
  }

  // get challenge
  pub fn fiat_shamir_prover (&self, num_bytes: usize) -> Bytes {
    let mut hasher = Shake256::default();

    let str = serde_json::to_string(&self.objects).unwrap();
    hasher.update(str.as_bytes());

    let mut buf = vec![0u8; num_bytes];
    let reader = hasher.finalize_xof();
    reader
      .take(num_bytes as u64)
      .read(&mut buf)
      .unwrap(); // under the hood it asserts so Result is useless

    Bytes::new(buf)
  }

  // reproduce challenge
  pub fn fiat_shamir_verifier (&self, num_bytes: usize) -> Bytes {
    let mut hasher = Shake256::default();

    let slice = &self.objects[0..self.read_index];
    let str = serde_json::to_string(slice).unwrap();
    hasher.update(str.as_bytes());

    let mut buf = vec![0u8; num_bytes];
    let reader = hasher.finalize_xof();
    reader
      .take(num_bytes as u64)
      .read(&mut buf)
      .unwrap(); // under the hood it asserts so Result is useless

    Bytes::new(buf)
  }
}

#[cfg(test)]
mod tests {
  use std::collections::HashMap;
  use super::*;

  #[test]
  fn order () {
    #[derive(Debug, Serialize, Clone)]
    enum SomeObjects {
      Map(HashMap<String, usize>),
      Vec(Vec<usize>),
      Str(String),
    }
    impl PartialEq for SomeObjects {
      fn eq (&self, other: &Self) -> bool {
        match self {
          SomeObjects::Map(m_1) => {
            match other {
              SomeObjects::Map(m_2) => m_1.eq(m_2),
              _ => false,
            }
          },
          SomeObjects::Vec(v_1) => {
            match other {
              SomeObjects::Vec(v_2) => v_1.eq(v_2),
              _ => false,
            }
          },
          SomeObjects::Str(s_1) => {
            match other {
              SomeObjects::Str(s_2) => s_1.eq(s_2),
              _ => false,
            }
          },
        }
      }
    }

    let mut proof_stream = ProofStream::<SomeObjects>::new();

    assert_eq!(proof_stream.fiat_shamir_prover(64), "ec784925b52067bce01fd820f554a34a3f8522b337f82e00ea03d3fa2b207ef9c2c1b9ed900cf2bbfcd19a232a94c6121e041615305c4155d46d52f58a8cff1c".into());

    proof_stream.push(SomeObjects::Str(String::from("Hello, World!")));
    proof_stream.push(SomeObjects::Vec(vec![0, 1, 5, 234]));
    let mut map = HashMap::new();
    map.insert(String::from("something"), 123);
    proof_stream.push(SomeObjects::Map(map.clone()));

    assert_eq!(proof_stream.fiat_shamir_prover(4), "5127a2a6".into());
    assert_eq!(proof_stream.fiat_shamir_prover(64), "5127a2a64ab1cbf4595c58f57171470e960d4d99f150be084a2410445916c8b325742d20aaa63f032915e64ac1c3556e93fb8d97d3abbb9d34f0a3636859b751".into());
    assert_eq!(proof_stream.fiat_shamir_verifier(64), "ec784925b52067bce01fd820f554a34a3f8522b337f82e00ea03d3fa2b207ef9c2c1b9ed900cf2bbfcd19a232a94c6121e041615305c4155d46d52f58a8cff1c".into());

    assert_eq!(proof_stream.pull(), Some(SomeObjects::Str(String::from("Hello, World!"))));
    assert_eq!(proof_stream.pull(), Some(SomeObjects::Vec(vec![0, 1, 5, 234])));
    assert_eq!(proof_stream.pull(), Some(SomeObjects::Map(map)));

    assert_eq!(proof_stream.fiat_shamir_verifier(64), "5127a2a64ab1cbf4595c58f57171470e960d4d99f150be084a2410445916c8b325742d20aaa63f032915e64ac1c3556e93fb8d97d3abbb9d34f0a3636859b751".into());
  }
}
