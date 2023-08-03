use std::ops::{BitAnd, Shr};
use std::usize;

#[derive(Debug)]
pub struct BitIter<T>(Option<usize>, T);

impl<T> BitIter<T> {
  pub fn bit_index (&self) -> Option<usize> {
    self.0
  }
}

impl Into<BitIter<u128>> for u128 {
  fn into (self) -> BitIter<u128> {
    if self != 0 {
      // for i in (0..u128::BITS).rev() {
      //   if (self >> i) & 1 == 1 {
      //     return BitIter(Some(i), self);
      //   }
      // }
      let pos = (0..u128::BITS)
        .rposition(|i| (self >> i) & 1 == 1);
      return BitIter(pos, self);
    }

    BitIter(None, 0)
  }
}

impl Into<BitIter<usize>> for usize {
  fn into (self) -> BitIter<usize> {
    if self != 0 {
      // for i in (0..usize::BITS).rev() {
      //   if (self >> i) & 1 == 1 {
      //     return BitIter(Some(i), self);
      //   }
      // }
      let pos = (0..usize::BITS)
        .rposition(|i| (self >> i) & 1 == 1);
      return BitIter(pos, self);
    }

    BitIter(None, 0)
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

impl<T> PartialEq for BitIter<T>
  where
    T: PartialEq
{
  fn eq (&self, other: &Self) -> bool {
    (self.1 == other.1) && (self.0 == other.0)
  }
}


#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_some () {
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
  fn test_none () {
    let mut s: BitIter<usize> = 0.into(); // 1011

    assert_eq!(s, BitIter(None, 0));
    assert_eq!(s.next(), None);
  }
}
