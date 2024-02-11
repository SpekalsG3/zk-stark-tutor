use std::ops::{BitXor, Rem, Shl};
use crate::field::field::Field;

#[derive(Debug, Clone, PartialEq)]
pub struct U512 {
  bits: [u128; 4]
}

impl U512 {
  pub fn new () -> Self {
    Self {
      bits: [0; 4],
    }
  }
}

impl From<u128> for U512 {
  fn from (value: u128) -> Self {
    Self {
      bits: [0, 0, 0, value],
    }
  }
}

impl Shl<u32> for U512 {
  type Output = Self;
  fn shl (mut self, rhs: u32) -> Self::Output {
    let mut overflows = [0_u128; 4];

    let full_overflow = rhs >= 128;
    let shift = rhs % 128;
    let overflow_target = rhs as usize / 128;
    for i in (0..4).rev() {
      let bit_i = self.bits.get_mut(i).unwrap();

      if overflow_target + 1 <= i {
        let target_i = i - overflow_target - 1;

        if overflows[target_i] == 0 {
          overflows[target_i] = bit_i.checked_shr(128 - shift).unwrap_or(0);
        }
      }

      if full_overflow {
        if overflow_target <= i {
          let target_i = i - overflow_target;

          if overflows[target_i] == 0 {
            overflows[target_i] = bit_i.checked_shl(shift).unwrap_or(0);
          }
        }

        *bit_i = overflows[i];
      } else {
        *bit_i = overflows[i] + bit_i.checked_shl(shift).unwrap_or(0);
      }
    }

    self
  }
}

impl BitXor<u128> for U512 {
  type Output = Self;
  fn bitxor (mut self, rhs: u128) -> Self::Output {
    let last = self.bits.last_mut().unwrap();
    *last = *last ^ rhs;
    self
  }
}

impl Rem<u128> for U512 {
  type Output = u128;
  fn rem (self, rhs: u128) -> Self::Output {
    let field = Field::new(rhs);

    self
      .bits
      .iter()
      .rev()
      .enumerate()
      .map(|(bit_offset, bit)| {
        let rem_res = bit % rhs;
        let rem_off = (0..bit_offset*2)
          .map(|_| (1 << (u128::BITS / 2)) % rhs)
          .fold(1, |acc, b| field.mul_mod(acc, b)); // (acc * b) % rhs
        field.mul_mod(rem_res, rem_off) // (rem_res * rem_off) % rhs
      })
      .fold(0, |acc, rem| field.add_mod(acc, rem)) // (acc + rem) % rhs
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn shl_less_128 () {
    let var: U512 = u128::MAX.into();
    assert_eq!(var.clone() << 1, U512 {
      bits: [0, 0, 1, u128::MAX << 1],
    });
    assert_eq!(var.clone() << 2, U512 {
      bits: [0, 0, 3, u128::MAX << 2],
    });

    let var = U512 {
      bits: [0b1, 0b0, u128::MAX, 0b1011],
    };
    assert_eq!(var << 4, U512 {
      bits: [0b10000, 0b1111, u128::MAX << 4, 0b10110000]
    });

    let var = U512 {
      bits: [0, 314322366231663367212220176324606862154, 84432643017724152666804417580452839161, 258876115543403587319877968525142574610]
    };
    assert_eq!(var << 8, U512 {
      bits: [236, 159887161964344628971957785202058807871, 176967496518259884509330632394529503682, 257506396449256441994086100673466077696],
    });
  }

  #[test]
  fn shl_more_128 () {
    let var: U512 = u128::MAX.into();
    assert_eq!(var.clone() << u128::BITS, U512 {
      bits: [0, 0, u128::MAX, 0],
    });
    assert_eq!(var.clone() << (u128::BITS + 1), U512 {
      bits: [0, 1, u128::MAX << 1, 0],
    });
    assert_eq!(var.clone() << u128::BITS * 2, U512 {
      bits: [0, u128::MAX, 0, 0],
    });
    assert_eq!(var.clone() << u128::BITS * 4, U512 {
      bits: [0, 0, 0, 0],
    });

    let var = U512 {
      bits: [0b10, 0b0, u128::MAX, 0b1011],
    };
    assert_eq!(var << (u128::BITS + 1), U512 {
      bits: [0b1, u128::MAX << 1, 0b10110, 0b0],
    });
  }

  #[test]
  fn xor () {
    let var = U512 {
      bits: [0, 0, 1, 0b1010],
    };
    assert_eq!(var ^ 0b0110, U512 {
      bits: [0, 0, 1, 0b1100],
    });

    let var = U512 {
      bits: [236, 159887161964344628971957785202058807871, 176967496518259884509330632394529503682, 257506396449256441994086100673466077696],
    };
    assert_eq!(var ^ 0b01010110, U512 {
      bits: [236, 159887161964344628971957785202058807871, 176967496518259884509330632394529503682, 257506396449256441994086100673466077782],
    })
  }

  #[test]
  fn rem () {
    let var = U512 {
      bits: [0, 0, 1, 0],
    };
    assert_eq!(var % 9, 4);

    let var = U512 {
      bits: [0, 1, 0, 0],
    };
    assert_eq!(var % 9, 7);

    let var = U512 {
      bits: [64587, 3827465, 13486, 1346],
    };
    assert_eq!(var % 18, 8);
  }
}
