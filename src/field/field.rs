use std::cmp::Ordering;
use serde::Serialize;
use crate::field::field_element::FieldElement;
use crate::utils::bytes::Bytes;
use crate::utils::u512::U512;

// 270497897142230380135924736767050121217
pub const FIELD_PRIME: u128 = 1 + 407 * (1 << 119);

#[derive(Debug, Serialize, Clone, PartialEq, PartialOrd)]
pub struct Field {
  pub order: u128,
}


impl Field {
  pub fn new (order: u128) -> Field {
    // assert_eq!(order, FIELD_PRIME, "Only 1+407*2^119 currently implemented");
    Field {
      order,
    }
  }
}

impl<'a> Field {
  pub fn generator (&'a self) -> FieldElement<'a> {
    assert_eq!(self.order, FIELD_PRIME, "Do not know generator for orders other than 1+407*2^119");
    FieldElement {
      field: &self,
      value: 85408008396924667383611388730472331217,
    }
  }

  pub fn primitive_nth_root (&'a self, n: u128) -> FieldElement<'a> {
    assert!((n & (n - 1) == 0) && (n <= (1 << 119)), "Field does not have any roots where n > 2^119 or not a power of two.");

    // same as generator, is it important?
    let mut root = self.generator();

    let mut order = 1 << 119;
    while order != n {
      root = root ^ 2;
      order = order / 2;
    }

    root
  }

  pub fn zero (&'a self) -> FieldElement<'a> {
    FieldElement {
      field: &self,
      value: 0,
    }
  }

  pub fn one (&'a self) -> FieldElement<'a> {
    FieldElement {
      field: &self,
      value: 1,
    }
  }

  pub fn sample (&'a self, bytes: &Bytes) -> FieldElement<'a> {
    let res = bytes
      .iter()
      .fold(U512::new(), |acc, b| {
        (acc << 8) ^ (*b as u128)
      });

    FieldElement {
      field: &self,
      value: res % self.order,
    }
  }

  pub(crate) fn sub_mod (&self, a: u128, b: u128) -> u128 {
    match a.cmp(&b) {
      Ordering::Greater => a - b,
      Ordering::Equal => 0,
      Ordering::Less => self.order - b + a,
    }
  }

  pub(crate) fn add_mod (&self, a: u128, b: u128) -> u128 {
    if b == 0 {
      return a;
    }

    self.sub_mod(a, self.order - b)
  }

  pub(crate) fn mul_mod (&self, a: u128, b: u128) -> u128 {
    let mut res = 0;

    let mut a = a;
    let mut b = b;

    while b > 0 {
      if b % 2 == 1 {
        res = self.add_mod(res, a);
      }

      a = self.add_mod(a, a);

      b /= 2;
    }

    // // source - https://www.youtube.com/watch?v=9hSmQtL49g4&ab_channel=DG
    // // doesnt work
    // while a != 0 {
    //   if a & 1 == 1  {
    //     res ^= b;
    //   }
    //
    //   a >>= 1;
    //   b <<= 1;
    //
    //   if degree(b) == self.degree {
    //     b ^= self.prime;
    //   }
    // }

    res
  }

  pub(crate) fn neg_mod (&self, a: u128) -> u128 {
    if a == 0 {
      0
    } else {
      self.order - a
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn mul () {
    let field = Field::new(FIELD_PRIME);

    assert_eq!(field.mul_mod(2, 3), 6);
    assert_eq!(field.mul_mod(FIELD_PRIME, 3), 0);
    assert_eq!(field.mul_mod(FIELD_PRIME - 1, 3), FIELD_PRIME - 3);
  }

  #[test]
  fn primitive_nth_root () {
    let field = Field::new(FIELD_PRIME);

    let n = 256;
    let n_log = 8; // so that 2 ** `n_log` = n
    let z = field.primitive_nth_root(n);

    assert_eq!(
      field.primitive_nth_root(256),
      FieldElement::new(&field, 178902808384765167578311106676137348214),
    );
    assert_eq!(
      field.primitive_nth_root(2),
      FieldElement::new(&field, 270497897142230380135924736767050121216),
    );

    // straightforward
    let powered = (0..n-1)
      .fold(z.value, |acc, _| {
        field.mul_mod(acc, z.value)
      });
    assert_eq!(powered, 1, "omega is not {}th root of unity", n);

    let powered = (0..n-2)
      .fold(z.value, |acc, _| {
        field.mul_mod(acc, z.value)
      });
    assert_ne!(powered, 1, "omega is not primitive");

    assert_eq!(z ^ (1 << n_log), field.one(), "omega not nth root of unity");
    assert_ne!(z ^ (1 << (n_log - 1)), field.one(), "omega not primitive");
  }

  #[test]
  fn sample () {
    let field = Field::new(FIELD_PRIME);

    let bytes: Bytes = "ec784925b52067bce01fd820f554a34a3f8522b337f82e00ea03d3fa2b207ef9c2c1b9ed900cf2bbfcd19a232a94c6121e041615305c4155d46d52f58a8cff1c"
      .into();
    assert_eq!(field.sample(&bytes), FieldElement {
      field: &field,
      value: 42271748005837835913754035451780029064,
    });

    let bytes = Bytes::from("6c9c4992");
    assert_eq!(field.sample(&bytes), FieldElement {
      field: &field,
      value: 1822181778,
    });

    let bytes = Bytes::from("ac4cd3be");
    assert_eq!(field.sample(&bytes), FieldElement {
      field: &field,
      value: 2890716094,
    });
  }

  #[test]
  fn neg () {
    let field = Field::new(FIELD_PRIME);
    assert_eq!(field.neg_mod(256), 270497897142230380135924736767050120961);
    let field = Field::new(100);
    assert_eq!(
      field.add_mod(20, field.neg_mod(20)),
      0,
    );
    assert_eq!(
      field.add_mod(20, field.neg_mod(19)),
      1,
    );
  }
}
