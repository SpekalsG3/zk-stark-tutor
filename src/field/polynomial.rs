use std::cmp::max;
use std::ops::{Add, BitXor, Div, Mul, Neg, Rem, Sub};
use crate::field::field_element::FieldElement;

#[derive(Debug, Clone)]
pub struct Polynomial<'a> {
  pub coefficients: Vec<FieldElement<'a>>,
}

impl<'a> Polynomial<'a> {
  pub fn new (coeffs: Vec<FieldElement<'a>>) -> Self {
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
        if c.eq(&zero) {
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
      None => None,
    }
  }

  pub fn evaluate (&self, point: &FieldElement<'a>) -> FieldElement<'a> {
    // todo check if coefficients are from the same field

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

  pub fn scale (self, factor: u128) -> Polynomial<'a> {
    // todo() need an example of calculation
    Polynomial::new(
      self
        .coefficients
        .iter()
        .enumerate()
        .map(|(i, el)| {
          // todo, i guess, has to be checked
          FieldElement {
            field: el.field,
            value: el.field.mul_mod(factor ^ i as u128, el.value),
          }
        })
        .collect()
    )
  }

  pub fn interpolate_domain <'m>(domain: &[FieldElement<'m>], values: &[FieldElement<'m>]) -> Polynomial<'m> {
    // todo() need an example of calculation
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
        ) * Polynomial::new(vec![ (*el_i - *el_j).inverse() ])
      }

      acc = acc + prod;
    }

    acc
  }

  pub fn zerofier_domain <'m>(domain: Vec<FieldElement<'m>>) -> Polynomial<'m> {
    // todo() need an example of calculation
    let field = domain.first().unwrap().field;
    let x = Polynomial::new(vec![field.zero(), field.one()]);

    domain
      .into_iter()
      .fold(Polynomial::new(vec![field.one()]), |acc, d| {
        acc * (x.clone() - Polynomial::new(vec![d]))
      })
  }

  pub fn test_colinearity <'m>(points: Vec<(FieldElement<'m>, FieldElement<'m>)>) -> bool {
    // todo() need an example of calculation
    let (domain, values) = points
      .into_iter()
      .fold((vec![], vec![]), |mut acc, p| {
        acc.0.push(p.0);
        acc.1.push(p.1);
        acc
      });

    let poly = Polynomial::interpolate_domain(&domain, &values);

    match poly.degree() {
      Some(d) => d <= 1,
      None => true,
    }
  }

  pub(crate) fn divide (numerator: Polynomial<'a>, denominator: Polynomial<'a>) -> Option<(Polynomial<'a>, Polynomial<'a>)> {
    let denom_degree = denominator.degree();
    if denom_degree.is_none() {
      return None;
    }
    let denom_degree = denom_degree.unwrap();
    let numer_degree = match numerator.degree() {
      Some(i) if i >= denom_degree => i,
      _ => {
        return Some((
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
    Some((quotient, remainder))
  }
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

impl<'a> PartialEq for Polynomial<'a> {
  fn eq (&self, other: &Self) -> bool {
    let degree = self.degree();
    if degree != other.degree() {
      return false;
    }
    if degree.is_none() {
      return true;
    }

    self
      .coefficients
      .iter()
      .enumerate()
      .all(|(i, el)| {
        other.coefficients.get(i).unwrap() == el
      })
  }

  fn ne (&self, other: &Self) -> bool {
    !self.eq(other)
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

impl<'a> Div for Polynomial<'a> {
  type Output = Self;
  fn div (self, rhs: Self) -> Self::Output {
    let res = Polynomial::divide(self, rhs);
    assert!(res.is_some(), "Denominator is empty or zero");

    let (q, r) = res.unwrap();
    assert!(r.is_zero(), "Cannot perform true division because remained is not zero");

    q
  }
}

impl<'a> Rem for Polynomial<'a> {
  type Output = Self;
  fn rem (self, rhs: Self) -> Self::Output {
    let res = Polynomial::divide(self, rhs);
    assert!(res.is_some(), "Denominator is empty or zero");

    let (_, r) = res.unwrap();

    r
  }
}

impl<'a> BitXor<u128> for Polynomial<'a> {
  type Output = Self;

  fn bitxor (self, exponent: u128) -> Self::Output {
    // todo: test, need an example of calculation
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

    for i in (0..128).rev() {
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

    let poly_c = poly_a.add(poly_b);

    assert_eq!(poly_c.coefficients, vec![
      FieldElement {
        field: &field,
        value: 7,
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
    assert_eq!(Polynomial::divide(nomin, denom), Some((quotient, remainder)));
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
        .map(|i| FieldElement {
          field: &field,
          value: i,
        })
        .collect::<Vec<_>>(),
      &vec![
        FieldElement {
          field: &field,
          value: 1,
        },
        FieldElement {
          field: &field,
          value: 4,
        },
        FieldElement {
          field: &field,
          value: 9,
        },
      ]
    );
    assert_eq!(poly, Polynomial::new(vec![
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
        value: 1,
      },
    ]));

    let domain = (1..7)
      .map(|i| FieldElement {
        field: &field,
        value: i,
      })
      .collect::<Vec<_>>();
    let values = vec![
      FieldElement {
        field: &field,
        value: 5,
      },
      FieldElement {
        field: &field,
        value: 2,
      },
      FieldElement {
        field: &field,
        value: 2,
      },
      FieldElement {
        field: &field,
        value: 1,
      },
      FieldElement {
        field: &field,
        value: 5,
      },
      FieldElement {
        field: &field,
        value: 0,
      },
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
      })
  }
}
