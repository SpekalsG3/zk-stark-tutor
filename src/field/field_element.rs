use std::fmt::{Debug, Formatter};
use std::ops::{Add, BitXor, Div, Mul, Neg, Sub};
use crate::field::field::Field;
use crate::utils::bit_iter::BitIter;
use crate::utils::bytes::Bytes;

#[derive(Clone, Copy, PartialEq, PartialOrd)]
pub struct FieldElement<'a> {
  pub field: &'a Field,
  pub value: u128,
}

impl Debug for FieldElement<'_> {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}", self.value)
  }
}

impl<'a> FieldElement<'a> {
  pub fn unstringify(str: &str, field: &'a Field) -> Self {
    Self {
      field,
      value: str.parse().unwrap(),
    }
  }

  pub fn new(field: &'a Field, value: u128) -> Self {
    Self {
      field,
      value,
    }
  }

  pub fn inverse (&self) -> FieldElement<'a> {
    FieldElement {
      field: self.field,
      value: self.field.inv(self.value),
    }
  }

  pub fn is_zero (self) -> bool {
    self.value == 0
  }
}

impl<'a> Into<Bytes> for FieldElement<'a> {
  fn into(self) -> Bytes {
    self.to_string().as_bytes().into()
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

    self * rhs.inverse()
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

// todo tutorial uses `BitXor` for power, replaces later
impl<'a> BitXor<usize> for FieldElement<'a> {
  type Output = Self;

  fn bitxor (self, exponent: usize) -> Self::Output {
    let mut acc = self.field.one();

    let iter: BitIter<usize> = exponent.into();
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
  fn mul () {
    let field = Field::new(FIELD_PRIME);
    assert_eq!(
      FieldElement::new(&field, 49789714223038013592473676705012096123) * FieldElement::new(&field, 6534789852937546098347957826345234),
      FieldElement::new(&field, 105250150227149389100670877502232671566),
    );

    let el_1 = FieldElement::new(&field, 8);
    let el_3 = FieldElement::new(&field, 12);
    assert_eq!(el_1 * el_3, FieldElement::new(&field, 96));

    let el_1 = FieldElement::new(&field, 3);
    let el_3 = FieldElement::new(&field, 270497897142230380135924736767050121215);
    assert_eq!(el_1 * el_3, FieldElement::new(&field, 270497897142230380135924736767050121211));
  }

  #[test]
  fn div () {
    let field = Field::new(FIELD_PRIME);
    assert_eq!(
      FieldElement::new(&field, 74658620945386735627456854792784352353) / FieldElement::new(&field, 85408008396924667383611388730472331217),
      FieldElement::new(&field, 120557879365253444230411244907275635216)
    );

    let el_1 = FieldElement::new(&field, 12);
    let el_2 = FieldElement::new(&field, 4);
    assert_eq!(el_1 / el_2, FieldElement::new(&field, 3));

    let el_1 = FieldElement::new(&field, 270497897142230380135924736767050121215);
    let el_2 = FieldElement::new(&field, 5);
    assert_eq!(el_1 / el_2, FieldElement::new(&field, 54099579428446076027184947353410024243));

    let el_1 = FieldElement::new(&field, 270497897142230380135924736767050121215);
    let el_2 = FieldElement::new(&field, 5);
    assert_eq!(el_1 / el_2, FieldElement::new(&field, 54099579428446076027184947353410024243));

    assert_eq!(
      FieldElement::new(&field, 5012096123) / FieldElement::new(&field, 6534789852937546098347957826345234),
      FieldElement::new(&field, 109071144973379706934869779239844248849),
    );

    let field = Field::new(8);
    let el_1 = FieldElement::new(&field, 2);
    let el_2 = FieldElement::new(&field, 7);
    assert_eq!(el_1 / el_2, FieldElement::new(
      &field,
      6, // because 6 * 7 = 2 (mod 8)
    ));
  }

  #[test]
  fn inverse () {
    let field = Field::new(FIELD_PRIME);

    assert_eq!(
      FieldElement::new(&field, 256).inverse(),
      FieldElement::new(&field, 269441264731518542713518780764053831681)
    );

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
    let field = Field::new(FIELD_PRIME);
    assert_eq!(
      FieldElement::new(&field, 270497897142230380135924736767050120961) + FieldElement::new(&field, 300),
      FieldElement::new(&field, 44),
    );

    let field = Field::new(100);
    assert_eq!(
      FieldElement::new(&field, 20) + FieldElement::new(&field, 20),
      FieldElement::new(&field, 40),
    );
    assert_eq!(
      FieldElement::new(&field, 20) + FieldElement::new(&field, 19).neg(),
      field.one(),
    );
    assert_eq!(
      FieldElement::new(&field, 80) + FieldElement::new(&field, 21),
      field.one(),
    );
  }

  #[test]
  fn sub () {
    let field = Field::new(FIELD_PRIME);
    assert_eq!(
      FieldElement::new(&field, 44) - FieldElement::new(&field, 200),
      FieldElement::new(&field, 270497897142230380135924736767050121061),
    );

    let field = Field::new(100);
    assert_eq!(
      FieldElement::new(&field, 20) - FieldElement::new(&field, 20),
      field.zero(),
    );
    assert_eq!(
      FieldElement::new(&field, 20) - FieldElement::new(&field, 19),
      field.one(),
    );
    assert_eq!(
      FieldElement::new(&field, 20) - FieldElement::new(&field, 21),
      field.one().neg(),
    );
  }

  #[test]
  fn neg () {
    let field = Field::new(FIELD_PRIME);
    assert_eq!(
      FieldElement::new(&field, 6534789852937546098).neg(),
      FieldElement::new(&field, 270497897142230380129389946914112575119),
    );

    let field = Field::new(100);
    assert_eq!(
      FieldElement::new(&field, 1).neg(),
      FieldElement::new(&field, 99),
    );
    assert_eq!(
      FieldElement::new(&field, 20).neg(),
      FieldElement::new(&field, 80),
    );
  }

  #[test]
  fn pow () {
    let field = Field::new(FIELD_PRIME);

    assert_eq!(FieldElement::new(&field, 6534789852937546098) ^ 501209126122_usize, FieldElement::new(&field, 256557788041265930815463337858691703671));

    let el = FieldElement::new(&field, 15);
    assert_eq!(el ^ 4_usize, FieldElement::new(&field, 50625));

    let el = FieldElement::new(&field, 270497897142230380135);
    assert_eq!(el ^ 8_usize, FieldElement::new(&field, 79016866124691016201920330826259043252));
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
