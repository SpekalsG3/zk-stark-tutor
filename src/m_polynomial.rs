use std::collections::HashMap;
use std::ops::{Add, BitXor, Mul, Neg, Sub};
use crate::field::field::Field;
use crate::field::field_element::FieldElement;
use crate::field::polynomial::Polynomial;
use crate::utils::bit_iter::BitIter;

type MPolynomialKey<'a> = HashMap<Vec<u128>, FieldElement<'a>>;

#[derive(Debug, Clone)]
pub struct MPolynomial<'a> {
  dictionary: MPolynomialKey<'a>,
}

impl<'a> MPolynomial<'a> {
  pub fn new (dict: MPolynomialKey<'a>) -> Self {
    // Multivariate polynomials are represented as dictionaries with exponent vectors
    // as keys and coefficients as values. E.g.:
    // f(x,y,z) = 17 + 2xy + 42z - 19 * x^6 * y^3 * z^12 is represented as:
    // {
    //     (0,0,0) => 17,
    //     (1,1,0) => 2,
    //     (0,0,1) => 42,
    //     (6,3,12) => -19,
    // }
    Self {
      dictionary: dict,
    }
  }

  pub fn zero () -> Self {
    Self {
      dictionary: HashMap::new(),
    }
  }

  pub fn constant (element: FieldElement<'a>) -> Self {
    let mut dict = MPolynomialKey::new();
    dict.insert(vec![0], element);

    MPolynomial {
      dictionary: dict
    }
  }

  // Returns the multivariate polynomials representing each indeterminates linear function
  // with a leading coefficient of one. For three indeterminates, returns:
  // [f(x,y,z) = x, f(x,y,z) = y, f(x,y,z) = z]
  pub fn variables (num_variables: usize, field: &'a Field) -> Vec<MPolynomial<'a>>{
    (0..num_variables)
      .map(|i| {
        let mut exponent = vec![0_u128; i];
        exponent.push(1);
        exponent.extend(vec![0_u128; num_variables - i - 1]);

        let mut dict = MPolynomialKey::new();
        dict.insert(exponent, field.one());

        MPolynomial::new(dict)
      })
      .collect()
  }

  pub fn lift (polynomial: Polynomial<'a>, variable_index: usize) -> MPolynomial<'a> {
    let mut acc = MPolynomial::zero();
    
    if polynomial.is_zero() {
      return acc;
    }
    
    let field = polynomial.coefficients.first().unwrap().field;
    let variables = MPolynomial::variables(variable_index + 1, &field);
    let x = variables.last().unwrap();
    
    for (i, el) in polynomial.coefficients.iter().enumerate() {
      acc = acc + MPolynomial::constant(*el) * (x.clone().bitxor(i as u128))
    }
    
    acc
  }

  pub fn is_zero (&self) -> bool {
    if self.dictionary.is_empty() {
      return true;
    }

    !self
      .dictionary
      .values()
      .find(|v| !v.is_zero())
      .is_some()
  }

  pub fn evaluate (self, point: Vec<FieldElement<'a>>) -> FieldElement<'a> {
    let mut acc = point.first().unwrap().field.zero();

    for (k, v) in self.dictionary {
      let mut prod = v;

      for (i, k) in k.iter().enumerate() {
        prod = prod * (point.get(i).unwrap().bitxor(*k));
      }

      acc = acc + prod;
    }

    acc
  }

  pub fn evaluate_symbolic (self, point: Vec<Polynomial<'a>>) -> Polynomial<'a> {
    let mut acc = Polynomial::zero();

    for (k, v) in self.dictionary {
      let mut prod = Polynomial::new(vec![v]);

      for (i, k) in k.iter().enumerate() {
        prod = prod * (point.get(i).unwrap().clone().bitxor(*k));
      }

      acc = acc + prod;
    }

    acc
  }
}

impl<'a> Neg for MPolynomial<'a> {
  type Output = Self;
  fn neg (self) -> Self::Output {
    let mut dict = MPolynomialKey::new();
    for (k, v) in self.dictionary {
      dict.insert(k, -v);
    }

    MPolynomial {
      dictionary: dict
    }
  }
}

// TODO: it's all untested, idk how
impl<'a> Add for MPolynomial<'a> {
  type Output = Self;
  fn add (self, rhs: Self) -> Self::Output {
    if self.dictionary.is_empty() {
      return rhs.clone();
    }
    if rhs.dictionary.is_empty() {
      return self.clone();
    }
    let mut dict = HashMap::new();

    let max_len_self = self
      .dictionary
      .keys()
      .max_by(|a, b| a.len().cmp(&b.len()))
      .unwrap()
      .len();
    let max_len_rhs = rhs
      .dictionary
      .keys()
      .max_by(|a, b| a.len().cmp(&b.len()))
      .unwrap()
      .len();
    let num_variables = max_len_self.max(max_len_rhs);

    for (k, v) in self.dictionary {
      let mut k = k;
      k.extend(vec![0; num_variables - k.len()]);

      dict.insert(k, v);
    }
    for (k, v) in rhs.dictionary {
      let mut k = k;
      k.extend(vec![0; num_variables - k.len()]);

      dict
        .entry(k)
        .and_modify(|dk| *dk = *dk + v)
        .or_insert(v);
    }

    MPolynomial {
      dictionary: dict,
    }
  }
}

impl<'a> Sub for MPolynomial<'a> {
  type Output = Self;
  fn sub (self, rhs: Self) -> Self::Output {
    self.add(-rhs)
  }
}

impl<'a> Mul for MPolynomial<'a> {
  type Output = Self;
  fn mul (self, rhs: Self) -> Self::Output {
    let max_len_self = self
      .dictionary
      .keys()
      .max_by(|a, b| a.len().cmp(&b.len()))
      .unwrap()
      .len();
    let max_len_rhs = rhs
      .dictionary
      .keys()
      .max_by(|a, b| a.len().cmp(&b.len()))
      .unwrap()
      .len();
    let num_variables = max_len_self.max(max_len_rhs);

    let mut dict = MPolynomialKey::new();
    for (k0, v0) in self.dictionary {
      for (k1, v1) in &rhs.dictionary {
        let mut exponent = vec![0; num_variables];

        for (i, v) in k0.iter().enumerate() {
          *(exponent.get_mut(i).unwrap()) += v;
        }
        for (i, v) in k1.iter().enumerate() {
          *(exponent.get_mut(i).unwrap()) += v;
        }

        dict
          .entry(exponent)
          .and_modify(|dk| *dk = *dk + v0 * v1.clone())
          .or_insert(v0 * v1.clone());
      }
    }

    MPolynomial::new(dict)
  }
}

impl<'a> BitXor<u128> for MPolynomial<'a> {
  type Output = Self;
  fn bitxor (self, rhs: u128) -> Self::Output {
    if self.is_zero() {
      return MPolynomial::zero();
    }

    let field = self
      .dictionary
      .values()
      .collect::<Vec<_>>()
      .first()
      .unwrap()
      .field;
    let num_variables = self
      .dictionary
      .keys()
      .collect::<Vec<_>>()
      .first()
      .unwrap()
      .len();
    let exp = vec![0; num_variables];

    let mut dict = MPolynomialKey::new();
    dict.insert(exp, field.one());

    let iter: BitIter<u128> = rhs.into();
    iter.fold(MPolynomial::new(dict), |acc, b| {
      let mut acc = acc.clone() * acc;
      if b {
        acc = acc * self.clone()
      }
      acc
    })
  }
}

#[cfg(test)]
mod tests {
  use crate::field::field::{Field, FIELD_PRIME};
  use super::*;
  
  #[test]
  fn test_is_zero () {
    let field = Field::new(FIELD_PRIME);

    let m_poly_constant = MPolynomial::constant(FieldElement {
      field: &field,
      value: 0,
    });
    assert!(m_poly_constant.is_zero());

    let m_poly_constant = MPolynomial::constant(FieldElement {
      field: &field,
      value: 1,
    });
    assert!(!m_poly_constant.is_zero());
  }
}

