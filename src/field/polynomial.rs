use std::cmp::max;
use std::fmt::{Debug, Formatter};
use std::ops::{Add, BitXor, Div, Mul, Neg, Rem, Sub};
use crate::field::field_element::FieldElement;
use crate::utils::bit_iter::BitIter;

#[derive(Clone, PartialEq, PartialOrd)]
pub struct Polynomial<'a> {
  pub coefficients: Vec<FieldElement<'a>>,
}

impl Debug for Polynomial<'_> {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    f.write_str("[")?;
    let mut iter = self.coefficients.iter();
    if let Some(c) = iter.next() {
      write!(f, "{:?}", c)?;
    }
    while let Some(c) = iter.next() {
      write!(f, ", {:?}", c)?;
    }
    f.write_str("]")
  }
}

impl<'a> Polynomial<'a> {
  pub fn new (coeffs: Vec<FieldElement<'a>>) -> Self {
    if coeffs.len() > 2 {
      let first = coeffs.first().unwrap();
      coeffs
        .iter()
        .skip(1)
        .for_each(|c| assert_eq!(first.field, c.field, "coefficients must be in the same field"))
    }

    Polynomial {
      coefficients: coeffs,
    }
  }

  pub fn zero () -> Self {
    Polynomial {
      coefficients: vec![],
    }
  }

  pub fn degree (&self) -> Option<usize> {
    if self.coefficients.len() == 0 {
      return None;
    }

    let zero = self.coefficients.first().unwrap().field.zero();
    self
      .coefficients
      .iter()
      .enumerate()
      .fold(None, |acc, (i, c)| {
        if c == &zero {
          acc
        } else {
          Some(i)
        }
      })
  }

  pub fn is_zero (&self) -> bool {
    self.degree().is_none()
  }

  pub fn leading_coefficient (&self) -> Option<&FieldElement<'a>> {
    match self.degree() {
      Some(i) => self.coefficients.get(i),
      None => self.coefficients.last(),
    }
  }

  pub fn evaluate (&self, point: &FieldElement<'a>) -> FieldElement<'a> {
    if self.coefficients.len() > 0 {
      assert_eq!(
        self.coefficients.first().unwrap().field,
        point.field,
        "cannot evaluate point in a different field",
      );
    }

    let (value, _) = self
      .coefficients
      .iter()
      .fold((
        point.field.zero(),
        point.field.one(),
      ), |(value, xi), c| {
        let res = (
          value + c.mul(xi),
          xi.mul(*point),
        );
        res
      });

    value
  }

  pub fn evaluate_domain (&self, domain: &[FieldElement<'a>]) -> Vec<FieldElement<'a>> {
    domain
      .iter()
      .map(|el| self.evaluate(el))
      .collect()
  }

  pub fn scale (&self, factor: FieldElement) -> Polynomial<'a> {
    Polynomial::new(
      self
        .coefficients
        .iter()
        .enumerate()
        .map(|(i, coef)| {
          let pow = factor ^ i;
          FieldElement {
            field: coef.field,
            value: coef.field.mul_mod(pow.value, coef.value),
          }
        })
        .collect()
    )
  }

  pub fn interpolate_domain <'m>(domain: &[FieldElement<'m>], values: &[FieldElement<'m>]) -> Polynomial<'m> {
    assert_eq!(domain.len(), values.len(), "number of elements in domain does not match number of values");
    assert!(domain.len() > 0, "Cannot interpolate between zero points");

    let field = domain.first().unwrap().field;
    let x = Polynomial::new(vec![field.zero(), field.one()]);
    let mut acc = Polynomial::new(vec![]);

    for (i, el_i) in domain.iter().enumerate() {
      let mut prod = Polynomial::new(vec![*values.get(i).unwrap()]);

      for (j, el_j) in domain.iter().enumerate() {
        if i == j {
          continue
        }

        prod = prod * (
          x.clone() - Polynomial::new(vec![*el_j])
        ) * Polynomial::new(vec![ (*el_i - *el_j).inverse() ]);
      }

      acc = acc + prod;
    }

    acc
  }

  pub fn zerofier_domain <'m>(domain: &[FieldElement<'m>]) -> Polynomial<'m> {
    let field = domain.first().unwrap().field;
    let x = Polynomial::new(vec![field.zero(), field.one()]);

    domain
      .into_iter()
      .fold(Polynomial::new(vec![field.one()]), |acc, d| {
        acc * (x.clone() - Polynomial::new(vec![*d]))
      })
  }

  pub fn test_colinearity <'m>(points: Vec<(FieldElement<'m>, FieldElement<'m>)>) -> bool {
    // todo need an example of calculation
    let (domain, values) = points
      .into_iter()
      .fold((vec![], vec![]), |mut acc, p| {
        acc.0.push(p.0);
        acc.1.push(p.1);
        acc
      });

    let poly = Polynomial::interpolate_domain(&domain, &values);

    match poly.degree() {
      Some(d) => d == 1,
      None => false,
    }
  }

  pub fn divide_with_rem(numerator: Polynomial<'a>, denominator: Polynomial<'a>) -> Result<(Polynomial<'a>, Polynomial<'a>), String> {
    let denom_degree = denominator.degree();
    if denom_degree.is_none() {
      return Err("Denominator is zero or empty".to_string());
    }
    let denom_degree = denom_degree.unwrap();
    let numer_degree = match numerator.degree() {
      Some(i) if i >= denom_degree => i,
      _ => {
        return Ok((
          Polynomial::zero(),
          numerator,
        ));
      }
    };

    let field = denominator.coefficients.first().unwrap().field;
    let mut remainder = numerator.clone();

    let steps = numer_degree - denom_degree + 1;
    let mut quotient_coefficients = vec![field.zero(); steps];

    let denom_lead = denominator.leading_coefficient().unwrap();
    for _ in 0..steps {
      let remainder_degree = match remainder.degree() {
        Some(i) if i >= denom_degree => i,
        _ => break, // None or (i < denom_degree)
      };

      let coefficient = remainder.leading_coefficient().unwrap().div(*denom_lead);
      let shift = remainder_degree - denom_degree;

      let mut subtrahend_coeffs = vec![field.zero(); shift];
      subtrahend_coeffs.push(coefficient);
      let subtrahend = Polynomial {
        coefficients: subtrahend_coeffs
      }.mul(denominator.clone());

      remainder = remainder.sub(subtrahend);

      *(quotient_coefficients.get_mut(shift).unwrap()) = coefficient;
    };

    let quotient = Polynomial {
      coefficients: quotient_coefficients,
    };
    Ok((quotient, remainder))
  }

// impl<'a> Div for Polynomial<'a> {
//   type Output = Self;
  pub fn div (self, rhs: Self) -> Result<Self, String> {
    let res = Polynomial::divide_with_rem(self, rhs)?;

    let (q, r) = res;
    if !r.is_zero() {
      return Err("Cannot perform true division because remained is not zero".to_string());
    }

    Ok(q)
  }
// }
}

impl<'a> Neg for Polynomial<'a> {
  type Output = Self;
  fn neg (self) -> Self::Output {
    Polynomial {
      coefficients: self
        .coefficients
        .into_iter()
        .map(|c| -c)
        .collect()
    }
  }
}

impl<'a> Add for Polynomial<'a> {
  type Output = Self;
  fn add (self, rhs: Self) -> Self::Output {
    if self.degree().is_none() {
      return rhs;
    }
    if rhs.degree().is_none() {
      return self;
    }

    let field = self.coefficients.first().unwrap().field;
    let mut coeffs = vec![field.zero(); max(self.coefficients.len(), rhs.coefficients.len())];

    for (i, c) in self.coefficients.into_iter().enumerate() {
      let res = coeffs.get(i).unwrap().add(c);
      *(coeffs.get_mut(i).unwrap()) = res;
    }

    for (i, c) in rhs.coefficients.into_iter().enumerate() {
      let res = coeffs.get(i).unwrap().add(c);
      *(coeffs.get_mut(i).unwrap()) = res;
    }

    Polynomial {
      coefficients: coeffs,
    }
  }
}

impl<'a> Sub for Polynomial<'a> {
  type Output = Self;
  fn sub (self, rhs: Self) -> Self::Output {
    self.add(-rhs)
  }
}

impl<'a> Mul for Polynomial<'a> {
  type Output = Self;
  fn mul (self, rhs: Self) -> Self::Output {
    if self.coefficients.len() == 0 || rhs.coefficients.len() == 0 {
      return Polynomial::zero();
    }

    let zero = self.coefficients.first().unwrap().field.zero();
    let mut buf = vec![zero; self.coefficients.len() + rhs.coefficients.len() - 1];

    for (i, self_el) in self.coefficients.iter().enumerate() {
      if self_el.is_zero() {
        continue // optimization for sparse polynomials
      }
      for (j, rhs_el) in rhs.coefficients.iter().enumerate() {
        let res = buf.get(i + j).unwrap().add(self_el.mul(*rhs_el));

        *(buf.get_mut(i + j).unwrap()) = res;
      }
    }

    Polynomial {
      coefficients: buf,
    }
  }
}

impl<'a> Rem for Polynomial<'a> {
  type Output = Self;
  fn rem (self, rhs: Self) -> Self::Output {
    let res = Polynomial::divide_with_rem(self, rhs).unwrap();

    let (_, r) = res;

    r
  }
}

// todo tutorial uses `BitXor` for power, replaces later
impl<'a> BitXor<u128> for Polynomial<'a> {
  type Output = Self;

  fn bitxor (self, exponent: u128) -> Self::Output {
    if self.is_zero() {
      return Polynomial::zero();
    }

    let one = self.coefficients.first().unwrap().field.one();
    let mut acc = Polynomial {
      coefficients: vec![one],
    };

    if exponent == 0 {
      return acc;
    }

    let iter: BitIter<u128> = exponent.into();
    for i in (0..iter.count()).rev() {
      acc = acc.clone().mul(acc);
      if (1 << i) & exponent != 0 {
        acc = acc.mul(self.clone());
      }
    }

    acc
  }
}

#[cfg(test)]
mod tests {
  use crate::field::field::{Field, FIELD_PRIME};
  use super::*;

  #[test]
  fn degree_none () {
    let field = Field::new(FIELD_PRIME);

    let poly = Polynomial {
      coefficients: vec![
        field.zero(),
        field.zero(),
      ],
    };

    assert_eq!(poly.degree(), None);
  }

  #[test]
  fn degree_one () {
    let field = Field::new(FIELD_PRIME);

    let poly = Polynomial {
      coefficients: vec![
        field.zero(),
        field.zero(),
        field.one(),
        field.zero(),
      ],
    };

    assert_eq!(poly.degree(), Some(2));
  }

  #[test]
  fn add () {
    let field = Field::new(FIELD_PRIME);

    let poly_a = Polynomial {
      coefficients: vec![
        FieldElement {
          field: &field,
          value: 5,
        },
        FieldElement {
          field: &field,
          value: 6,
        },
      ],
    };
    let poly_b = Polynomial {
      coefficients: vec![
        FieldElement {
          field: &field,
          value: 2,
        },
      ],
    };

    assert_eq!(poly_a.sub(poly_b).coefficients, vec![
      FieldElement {
        field: &field,
        value: 3,
      },
      FieldElement {
        field: &field,
        value: 6,
      },
    ]);
  }

  #[test]
  fn divide () {
    let field = Field::new(FIELD_PRIME);

    let nomin = Polynomial {
      coefficients: vec![
        FieldElement {
          field: &field,
          value: 0,
        },
        FieldElement {
          field: &field,
          value: 1,
        },
        FieldElement {
          field: &field,
          value: 3,
        },
        FieldElement {
          field: &field,
          value: 18,
        },
        FieldElement {
          field: &field,
          value: 6,
        },
      ]
    };
    let denom = Polynomial {
      coefficients: vec![
        FieldElement {
          field: &field,
          value: 0,
        },
        FieldElement {
          field: &field,
          value: 0,
        },
        FieldElement {
          field: &field,
          value: 3,
        },
      ]
    };
    let quotient = Polynomial {
      coefficients: vec![
        FieldElement {
          field: &field,
          value: 1,
        },
        FieldElement {
          field: &field,
          value: 6,
        },
        FieldElement {
          field: &field,
          value: 2,
        },
      ]
    };
    let remainder = Polynomial {
      coefficients: vec![
        FieldElement {
          field: &field,
          value: 0,
        },
        FieldElement {
          field: &field,
          value: 1,
        },
        FieldElement {
          field: &field,
          value: 0,
        },
        FieldElement {
          field: &field,
          value: 0,
        },
        FieldElement {
          field: &field,
          value: 0,
        },
      ]
    };
    assert_eq!(Polynomial::divide_with_rem(nomin, denom), Ok((quotient, remainder)));
  }

  #[test]
  fn evaluate () {
    let field = Field::new(FIELD_PRIME);

    let poly = Polynomial {
      coefficients: vec![
        FieldElement {
          field: &field,
          value: 5,
        },
        FieldElement {
          field: &field,
          value: 0,
        },
        FieldElement {
          field: &field,
          value: 10,
        },
      ]
    };
    let point = FieldElement {
      field: &field,
      value: 3,
    };
    let res = FieldElement {
      field: &field,
      value: 95,
    };
    assert_eq!(poly.evaluate(&point), res);
  }

  #[test]
  fn interpolate () {
    let field = Field::new(FIELD_PRIME);

    let poly = Polynomial::interpolate_domain(
      &(1..4)
        .map(|i| FieldElement::new(&field, i))
        .collect::<Vec<_>>(),
      &vec![
        FieldElement::new(&field, 1),
        FieldElement::new(&field, 4),
        FieldElement::new(&field, 9),
      ]
    );
    assert_eq!(poly, Polynomial::new(vec![
      FieldElement::new(&field, 0),
      FieldElement::new(&field, 0),
      FieldElement::new(&field, 1),
    ]));

    let domain = (1..7)
      .map(|i| FieldElement::new(&field, i))
      .collect::<Vec<_>>();
    let values = vec![
      FieldElement::new(&field, 5),
      FieldElement::new(&field, 2),
      FieldElement::new(&field, 2),
      FieldElement::new(&field, 1),
      FieldElement::new(&field, 5),
      FieldElement::new(&field, 0),
    ];
    let poly = Polynomial::interpolate_domain(
      &domain,
      &values,
    );
    values
      .iter()
      .enumerate()
      .for_each(|(i, v)| {
        assert_eq!(&poly.evaluate(domain.get(i).unwrap()), v)
      });
    assert_ne!(poly.evaluate(&FieldElement {
      field: &field,
      value: 363,
    }), field.zero());
    assert_eq!(poly.degree(), Some(domain.len() - 1));
  }

  #[test]
  fn pow () {
    let field = Field::new(FIELD_PRIME);

    let poly = Polynomial {
      coefficients: vec![
        FieldElement {
          field: &field,
          value: 2,
        },
        FieldElement {
          field: &field,
          value: 5,
        },
      ],
    };

    assert_eq!(poly ^ 2, Polynomial {
      coefficients: vec![
        FieldElement {
          field: &field,
          value: 4,
        },
        FieldElement {
          field: &field,
          value: 20,
        },
        FieldElement {
          field: &field,
          value: 25,
        },
      ]
    })
  }

  #[test]
  fn scale () {
    let field = Field::new(FIELD_PRIME);

    let poly = Polynomial::new(vec![
      FieldElement::new(&field, 10),
      FieldElement::new(&field, 345),
      FieldElement::new(&field, 0),
      FieldElement::new(&field, 65),
      FieldElement::new(&field, 74),
      FieldElement::new(&field, 5),
    ]);
    let poly_res = Polynomial::new(vec![
      FieldElement::new(&field, 10),
      FieldElement::new(&field, 1380),
      FieldElement::new(&field, 0),
      FieldElement::new(&field, 4160),
      FieldElement::new(&field, 18944),
      FieldElement::new(&field, 5120),
    ]);
    assert_eq!(poly.scale(FieldElement::new(&field, 4)), poly_res);
  }

  #[test]
  fn zerofier () {
    let field = Field::new(FIELD_PRIME);

    let domain = vec![
      FieldElement::new(&field, 10),
      FieldElement::new(&field, 345),
      FieldElement::new(&field, 0),
      FieldElement::new(&field, 65),
      FieldElement::new(&field, 74),
      FieldElement::new(&field, 5),
    ];
    let zerofier = Polynomial::zerofier_domain(&domain);
    for d in domain {
      assert_eq!(zerofier.evaluate(&d), field.zero())
    }
  }

  // todo
  // #[test]
  // fn neg () {
  //   unimplemented!()
  // }

  // todo
  // #[test]
  // fn sub () {
  //   unimplemented!()
  // }

  // todo
  // #[test]
  // fn mul () {
  //   unimplemented!()
  // }

  // todo
  // #[test]
  // fn rem () {
  //   unimplemented!()
  // }
}
