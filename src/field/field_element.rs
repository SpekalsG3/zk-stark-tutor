use std::cmp::Ordering;
use std::ops::{Add, BitXor, Div, Mul, Neg, Sub};
use serde::Serialize;
use crate::field::field::Field;
use crate::utils::bit_iter::BitIter;
use crate::utils::xgcd::u_xgcd;

#[derive(Debug, Clone, Copy, Serialize, PartialEq, PartialOrd)]
pub struct FieldElement<'a> {
  pub field: &'a Field,
  pub value: u128,
}

impl<'a> FieldElement<'a> {
  // inverse of `x` is `x ** -1 = 1/x` so that `x` multiplied by inversed `x` is `1`
  pub fn inverse (&self) -> FieldElement<'a> {
    let (a, _, _) = u_xgcd(self.value, self.field.order);

    // because a can be negative
    let a = match a.cmp(&0) {
      Ordering::Greater => a as u128,
      Ordering::Equal => 0,
      Ordering::Less => self.field.sub_mod(self.field.order, a.neg() as u128),
    };
    
    FieldElement {
      field: self.field,
      value: a,
    }
  }

  pub fn bytes (&self) -> Vec<u8> {
    self.to_string().into_bytes()
  }

  pub fn is_zero (self) -> bool {
    self.value == 0
  }
}

impl<'a> Add for FieldElement<'a> {
  type Output = Self;
  fn add (self, rhs: Self) -> Self::Output {
    FieldElement {
      field: self.field,
      value: self.field.add_mod(self.value, rhs.value),
    }
  }
}

impl<'a> Sub for FieldElement<'a> {
  type Output = Self;
  fn sub (self, rhs: Self) -> Self::Output {
    FieldElement {
      field: self.field,
      value: self.field.sub_mod(self.value, rhs.value),
    }
  }
}

impl<'a> Mul for FieldElement<'a> {
  type Output = Self;
  fn mul (self, rhs: Self) -> Self::Output {
    FieldElement {
      field: self.field,
      value: self.field.mul_mod(self.value, rhs.value),
    }
  }
}

impl<'a> Div for FieldElement<'a> {
  type Output = FieldElement<'a>;
  fn div (self, rhs: Self) -> Self::Output {
    assert_ne!(rhs.value, 0, "divide by zero");

    // assertion for local usage, behaviour that I don't understand
    let (a, _, _) = u_xgcd(rhs.value, self.field.order);

    // to prevent loosing data during cast
    let value: u128 = match a.cmp(&0) {
      Ordering::Greater => self.field.mul_mod(self.value, a as u128),
      Ordering::Less => self.field.mul_mod(self.value, self.field.neg_mod((-a) as u128)),
      Ordering::Equal => 0,
    };

    // assert_eq!(self.value / rhs.value, value, "division with remainder is not consistent");

    FieldElement {
      field: self.field,
      value: value,
    }
  }
}

impl<'a> Neg for FieldElement<'a> {
  type Output = FieldElement<'a>;
  fn neg (self) -> Self::Output {
    FieldElement {
      field: self.field,
      value: self.field.neg_mod(self.value),
    }
  }
}

impl<'a> ToString for FieldElement<'a> {
  fn to_string (&self) -> String {
    self.value.to_string()
  }
}

// todo tutorial uses `BitXor` for power, replaces later
impl<'a> BitXor<u128> for FieldElement<'a> {
  type Output = Self;

  fn bitxor (self, exponent: u128) -> Self::Output {
    let mut acc = self.field.one();

    let iter: BitIter<u128> = exponent.into();
    for i in (0..iter.count()).rev() {
      acc = acc * acc;
      if ((1 << i) & exponent) != 0 {
        acc = acc * self;
      }
    }
    
    acc
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::field::field::FIELD_PRIME;

  #[test]
  fn test_mul () {
    let field = Field::new(FIELD_PRIME);

    let el_1 = FieldElement {
      field: &field,
      value: 8,
    };
    let el_3 = FieldElement {
      field: &field,
      value: 12,
    };
    assert_eq!(el_1 * el_3, FieldElement {
      field: &field,
      value: 96,
    });

    let el_1 = FieldElement {
      field: &field,
      value: 3,
    };
    let el_3 = FieldElement {
      field: &field,
      value: 270497897142230380135924736767050121215,
    };
    assert_eq!(el_1 * el_3, FieldElement {
      field: &field,
      value: 270497897142230380135924736767050121211,
    });
  }

  #[test]
  fn test_div_without_rem () {
    let field = Field::new(FIELD_PRIME);

    let el_1 = FieldElement {
      field: &field,
      value: 12,
    };
    let el_2 = FieldElement {
      field: &field,
      value: 4,
    };
    assert_eq!(el_1 / el_2, FieldElement {
      field: &field,
      value: 3,
    });

    let el_1 = FieldElement {
      field: &field,
      value: 270497897142230380135924736767050121215,
    };
    let el_2 = FieldElement {
      field: &field,
      value: 5,
    };
    assert_eq!(el_1 / el_2, FieldElement {
      field: &field,
      value: 54099579428446076027184947353410024243,
    });

    let el_1 = FieldElement {
      field: &field,
      value: 270497897142230380135924736767050121215,
    };
    let el_2 = FieldElement {
      field: &field,
      value: 5,
    };
    assert_eq!(el_1 / el_2, FieldElement {
      field: &field,
      value: 54099579428446076027184947353410024243,
    });
  }

  #[test]
  fn test_div_with_rem () {
    let field = Field {
      order: 8,
    };
    let el_1 = FieldElement {
      field: &field,
      value: 2,
    };
    let el_2 = FieldElement {
      field: &field,
      value: 7,
    };
    assert_eq!(el_1 / el_2, FieldElement {
      field: &field,
      value: 6, // because 6 * 7 = 2 (mod 8) 
    });
  }

  #[test]
  fn test_inverse () {
    let field = Field::new(FIELD_PRIME);

    let el = FieldElement {
      field: &field,
      value: 8,
    };
    assert_eq!(el * el.inverse(), field.one());

    let el = FieldElement {
      field: &field,
      value: 270497897142230380135924736767050121215,
    };
    assert_eq!(el * el.inverse(), field.one());
  }

  #[test]
  fn add () {
    let field = Field::new(100);
    assert_eq!(
      FieldElement {
        field: &field,
        value: 20,
      } + FieldElement {
        field: &field,
        value: 20,
      },
      FieldElement {
        field: &field,
        value: 40,
      },
    );
    assert_eq!(
      FieldElement {
        field: &field,
        value: 20,
      } + (FieldElement {
        field: &field,
        value: 19,
      }).neg(),
      field.one(),
    );
    assert_eq!(
      FieldElement {
        field: &field,
        value: 80,
      } + FieldElement {
        field: &field,
        value: 21,
      },
      field.one(),
    );
  }

  #[test]
  fn sub () {
    let field = Field::new(100);
    assert_eq!(
      FieldElement {
        field: &field,
        value: 20,
      } - FieldElement {
        field: &field,
        value: 20,
      },
      field.zero(),
    );
    assert_eq!(
      FieldElement {
        field: &field,
        value: 20,
      } - FieldElement {
        field: &field,
        value: 19,
      },
      field.one(),
    );
    assert_eq!(
      FieldElement {
        field: &field,
        value: 20,
      } - FieldElement {
        field: &field,
        value: 21,
      },
      field.one().neg(),
    );
  }

  #[test]
  fn neg () {
    let field = Field::new(100);
    assert_eq!(FieldElement {
        field: &field,
        value: 1,
      }.neg(),
       FieldElement {
         field: &field,
         value: 99,
       },
    );
    assert_eq!(FieldElement { 
        field: &field, 
        value: 20, 
      }.neg(),
       FieldElement {
         field: &field,
         value: 80,
       },
    );
  }

  #[test]
  fn pow () {
    let field = Field::new(FIELD_PRIME);

    let el = FieldElement {
      field: &field,
      value: 15,
    };
    assert_eq!(el ^ 4, FieldElement {
      field: &field,
      value: 50625,
    });

    let el = FieldElement {
      field: &field,
      value: 270497897142230380135,
    };
    assert_eq!(el ^ 8, FieldElement {
      field: &field,
      value: 79016866124691016201920330826259043252,
    });
  }

  // fn bitxor () {
  //   let field = Field::new(FIELD_PRIME);
  // 
  //   let el = FieldElement {
  //     field: &field,
  //     value: 0b110011010, // 410
  //   };
  //   assert_eq!(el ^ 0b101101, FieldElement { // 45
  //     field: &field,
  //     value: 0b110110111, // 439
  //   });
  // }
}
