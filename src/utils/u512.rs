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
    // rhs = 1
    // i 0: bits_to_move = 0
    // i 1: bits_to_move = 0
    // i 2: bits_to_move = 0
    // i 3: bits_to_move = 1
    // rhs = 511
    // i 0: bits_to_move = 127 == rhs - ((3 - i) * 128)
    //      bits[0] - bits_to_move
    // i 1: bits_to_move = 128 == rhs - ((3 - i) * 128)
    // i 2: bits_to_move = 128
    // i 3: bits_to_move = 128

    // [0001, 1111, 1111, 1111] << 2
    // iter 1
    // move first chunk - [0100, 1111, 1111, 1111]
    // moved = 2
    // iter 2
    // subtract mask    - [0100, 0011, 1111, 1111]
    // add mask to prev - [0111, 0011, 1111, 1111]
    // move current     - [0111, 1100, 1111, 1111]

    // [0000, 0000, 0011, 1111] << 2
    // overflow = 0
    // iter 1  =>  i == 3
    // calc next_overflow = 1111 >> 2 = 11
    // set leftover  = 1111 << 2 = 1100  + overflow == 1100  =>  [0000, 0000, 0011, 1100]
    // set overflow = next_overflow
    // iter 2  => i == 2
    // calc next_overflow = 0011 >> 2 = 0
    // set leftover  = 0011 << 2 = 1100  + overflow == 1111  =>  [0000, 0000, 1111, 1100]
    // set overflow = next_overflow
    // rhs = 1   => shift = 1
    // rhs = 127 => shift = 127
    // rhs = 128 => shift = 128
    // rhs = 129 => shift = 128
    // rhs = 512 => shift = 128

    let mut overflows = [0_u128; 4];

    let (shift, overflow_target) = if rhs > 128 {
      panic!("not ready") // todo, save another half of `bit_i.checked_shr` to overflows[.0 + 1]
      // let offset = rhs / 128 + 1;
      // (128, (offset, Some(offset + 1)))
    } else {
      (rhs, (1, None::<usize>))
    };
    for i in (0..4).rev() {
      let mut bit_i = self.bits.get_mut(i).unwrap();

      if overflow_target.0 <= i {
        let target_i = i - overflow_target.0;
        overflows[target_i] = match bit_i.checked_shr(128 - shift) {
          Some(u) => u,
          None => 0,
        };

        match overflow_target.1 {
          Some(i) => {
            // overflows[overflow_target.1] = match bit_i.checked_shl(128 - shift) {
            //   Some(u) => u,
            //   None => 0,
            // };
          },
          None => {},
        }
      }

      *bit_i = overflows[i] + match bit_i.checked_shl(shift) {
        Some(u) => u,
        None => 0,
      };

      // {
      //   let mut bit_i = self.bits.get_mut(i).unwrap();
      //   let next_overflow = match bit_i.checked_shr(rhs) {
      //     Some(u) => u,
      //     None => 0,
      //   };
      //   println!("overflow: {}", overflow);
      //   *bit_i = overflow + match bit_i.checked_shl(rhs) {
      //     Some(u) => u,
      //     None => 0,
      //   };
      //   overflow = next_overflow
      // }

      // {
      //   let mut bit_prev_i = self.bits.get_mut(i - 1).unwrap();
      //   *bit_prev_i = *bit_prev_i + overflow;
      // }
      
      // let offset = (3 - i as u128) * 128 - last_bits;
      // if offset > rhs {
      //   continue;
      // }
      // 
      // let bits_to_move = rhs - offset;
      // let value = (0..bits_to_move)
      //   .fold(0, |acc, i| {
      //     acc + (1 << i)
      //   });
      // last_bits = bits_to_move;
      // println!("value: {}", value);
      // 
      // {
      //   let mut bit_i = self.bits.get_mut(i).unwrap();
      //   *bit_i = *bit_i << bits_to_move;
      // }
      // 
      // if i > 0 {
      //   let mut bit_prev_i = self.bits.get_mut(i - 1).unwrap();
      //   *bit_prev_i += value; 
      // }
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
  fn shl () {
    // let x: u128 = 340282366920938463463374607431768211455;
    // println!("x << 1 == {}", x << 64);
    // println!("x << 1 == {:?}", x.checked_shl(127));

    let var: U512 = u128::MAX.into();
    assert_eq!(var.clone() << 1, U512 {
      bits: [0, 0, 1, u128::MAX - 1],
    });
    assert_eq!(var.clone() << 2, U512 {
      bits: [0, 0, 3, u128::MAX - 3],
    });
    assert_eq!(var.clone() << u128::BITS, U512 {
      bits: [0, 0, u128::MAX, 0],
    });
    // assert_eq!(var.clone() << 256, U512 {
    //   bits: [0, u128::MAX, 0, 0],
    // });
    // assert_eq!(var.clone() << 512, U512 {
    //   bits: [0, 0, 0, 0],
    // });

    let var = U512 {
      bits: [0b1, 0b0, u128::MAX, 0b1011],
    };
    assert_eq!(var << 4, U512 {
      bits: [0b10000, 0b1111, u128::MAX - 15, 0b10110000]
    });

    let var = U512 {
      bits: [0, 314322366231663367212220176324606862154, 84432643017724152666804417580452839161, 258876115543403587319877968525142574610]
    };
    assert_eq!(var << 8, U512 {
      bits: [236, 159887161964344628971957785202058807871, 176967496518259884509330632394529503682, 257506396449256441994086100673466077696],
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
