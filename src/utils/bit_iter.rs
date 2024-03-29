use std::ops::{BitAnd, Shr};
use std::usize;

#[derive(Debug, PartialOrd, PartialEq)]
pub struct BitIter<T>(pub(crate) Option<usize>, pub(crate) T);

impl<T> BitIter<T> {
  pub fn bit_index (&self) -> Option<usize> {
    self.0
  }
}

impl<T> BitIter<T>
  where
    Self: Iterator
{
  pub fn degree (self) -> usize {
    if let Some((i, _)) = self.enumerate().last() {
      i
    } else {
      0
    }
  }
}

impl From<u128> for BitIter<u128> {
  fn from (value: u128) -> Self {
    if value == 0 {
      return BitIter(Some(0), 0);
    }

    // for i in (0..u128::BITS).rev() {
    //   if (value >> i) & 1 == 1 {
    //     return BitIter(Some(i), value);
    //   }
    // }
    let pos = (0..u128::BITS)
      .rposition(|i| (value >> i) & 1 == 1);
    return BitIter(pos, value);
  }
}

impl From<usize> for BitIter<usize> {
  fn from (value: usize) -> Self {
    if value == 0 {
      return BitIter(Some(0), 0);
    }

    // for i in (0..usize::BITS).rev() {
    //   if (value >> i) & 1 == 1 {
    //     return BitIter(Some(i), value);
    //   }
    // }
    let pos = (0..usize::BITS)
      .rposition(|i| (value >> i) & 1 == 1);
    BitIter(pos, value)
  }
}

impl<T> Iterator for BitIter<T>
  where
    T: Shr<usize> + Copy + From<u8>,
    <T as Shr<usize>>::Output: BitAnd<T>,
    <<T as Shr<usize>>::Output as BitAnd<T>>::Output: PartialEq<T>,
{
  type Item = bool;
  fn next (&mut self) -> Option<Self::Item> {
    let i = match self.0 {
      Some(i) => i,
      None => {
        return None;
      }
    };

    let bit = self.1 >> i & 1_u8.into();

    if i == 0 {
      self.0 = None;
    } else {
      self.0 = Some(i - 1);
    }

    Some(bit.eq(&1_u8.into()))
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn degree () {
    assert_eq!(BitIter::from(0b100_usize).degree(), 2);
    assert_eq!(BitIter::from(0b1010_usize).degree(), 3);
    assert_eq!(BitIter::from(0b1_u128 << 119).degree(), 119);
  }

  #[test]
  fn some () {
    let mut s: BitIter<usize> = 11.into(); // 1011

    assert_eq!(s, BitIter(Some(3), 11));
    assert_eq!(s.next(), Some(true));
    assert_eq!(s.next(), Some(false));
    assert_eq!(s.next(), Some(true));
    assert_eq!(s.next(), Some(true));
    assert_eq!(s.next(), None);
    assert_eq!(s, BitIter(None, 11));
  }

  #[test]
  fn zero () {
    let mut s: BitIter<usize> = 0.into(); // 0

    assert_eq!(s, BitIter(Some(0), 0));
    assert_eq!(s.next(), Some(false));
    assert_eq!(s.next(), None);
    assert_eq!(s, BitIter(None, 0));
  }
}
