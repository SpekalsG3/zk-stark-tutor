use std::collections::HashMap;
use std::ops::{Add, BitXor, Mul, Neg, Sub};
use crate::field::field::Field;
use crate::field::field_element::FieldElement;
use crate::field::polynomial::Polynomial;
use crate::utils::bit_iter::BitIter;

// todo: remake to `Vec<(FieldElement, Vec<u128>)>` - more efficient and less expensive
type MPolynomialKey<'a> = HashMap<Vec<u128>, FieldElement<'a>>;

#[derive(Debug, Clone)]
pub struct MPolynomial<'a> {
  pub(crate) dictionary: MPolynomialKey<'a>,
}

impl<'a> MPolynomial<'a> {
  pub fn new (dictionary: MPolynomialKey<'a>) -> Self {
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
      dictionary,
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
        let mut exponent = vec![0; i];
        exponent.push(1);
        exponent.resize(num_variables, 0);

        let mut dict = MPolynomialKey::new();
        dict.insert(exponent, field.one());

        MPolynomial::new(dict)
      })
      .collect()
  }

  pub fn lift (polynomial: &Polynomial<'a>, variable_index: usize) -> MPolynomial<'a> {
    let mut acc = MPolynomial::zero();

    if polynomial.is_zero() {
      return acc;
    }

    let field = polynomial.coefficients.first().unwrap().field;
    let variables = MPolynomial::variables(variable_index + 1, &field);
    let x = variables.last().unwrap();

    for (i, el) in polynomial.coefficients.iter().enumerate() {
      acc = acc + MPolynomial::constant(*el) * (x.clone() ^ i as u128)
    }

    acc
  }

  pub fn is_zero (&self) -> bool {
    if self.dictionary.is_empty() {
      return true;
    }

    !self.dictionary
      .values()
      .find(|v| !v.is_zero())
      .is_some()
  }

  pub fn evaluate (&self, point: &[FieldElement<'a>]) -> FieldElement<'a> {
    // self.dictionary
    //     .iter()
    //     .fold(acc, |acc, (exponents, coeff)| {
    //       let prod = exponents
    //           .iter()
    //           .enumerate()
    //           .fold(*coeff, |coeff, (index, exponent)| {
    //             let point_value = *point.get(index).unwrap();
    //             coeff * (point_value ^ *exponent)
    //           });
    //       acc + prod
    //     })

    let mut acc = point.first().unwrap().field.zero();

    for (exponents, coeff) in &self.dictionary {
      let mut prod = *coeff;

      for (index, exponent) in exponents.iter().enumerate() {
        let point_value = *point.get(index).unwrap();
        prod = prod * (point_value ^ *exponent);
      }

      acc = acc + prod;
    }

    acc
  }

  pub fn evaluate_symbolic (&self, point: &[Polynomial<'a>]) -> Polynomial<'a> {
    let mut acc = Polynomial::zero();

    for (exponents, coeff) in &self.dictionary {
      let mut prod = Polynomial::new(vec![*coeff]);

      for (index, exponent) in exponents.iter().enumerate() {
        let point_value = point.get(index).unwrap().clone();
        prod = prod * (point_value ^ *exponent);
      }

      acc = acc + prod;
    }

    acc
  }
}

impl<'a> PartialEq for MPolynomial<'a> {
  fn eq(&self, other: &Self) -> bool {
    if self.dictionary.len() != other.dictionary.len() {
      return false;
    }

    for (k_self, v_self) in &self.dictionary {
      match other.dictionary.get(k_self) {
        None => {
          return false;
        }
        Some(v_other) => {
          if v_self != v_other {
            return false
          }
        }
      }
    }

    return true
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

    let max_len_self = self.dictionary
      .keys()
      .max_by(|a, b| a.len().cmp(&b.len()))
      .unwrap()
      .len();
    let max_len_rhs = rhs.dictionary
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
    let max_len_self = self.dictionary
      .keys()
      .max_by(|a, b| a.len().cmp(&b.len()))
      .unwrap()
      .len();
    let max_len_rhs = rhs.dictionary
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

// todo: tutorial uses bitxor operator as power - replace with `pow` method
impl<'a> BitXor<u128> for MPolynomial<'a> {
  type Output = Self;
  fn bitxor (self, rhs: u128) -> Self::Output {
    if self.is_zero() {
      return MPolynomial::zero();
    }

    let field = self.dictionary
      .values()
      .next()
      .unwrap()
      .field;
    let num_variables = self.dictionary
      .keys()
      .next()
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
  fn is_zero () {
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

  #[test]
  fn mul () {
    let field = Field::new(FIELD_PRIME);

    {
      let mut dictionary = HashMap::new();
      dictionary.insert(vec![0, 1, 5], FieldElement::new(&field, 17));
      dictionary.insert(vec![42, 1, 5], FieldElement::new(&field, 5));
      // 17 * y * z^5 + 5 * x^42 * y * z^5
      let poly_a = MPolynomial { dictionary };

      let mut dictionary = HashMap::new();
      dictionary.insert(vec![42, 0], FieldElement::new(&field, 8));
      dictionary.insert(vec![0, 0], FieldElement::new(&field, field.neg_mod(7)));
      // 8 * x^42 - 7
      let poly_b = MPolynomial { dictionary };

      let mut dictionary = HashMap::new();
      dictionary.insert(vec![42, 1, 5], FieldElement::new(&field, 136) + FieldElement::new(&field, field.mul_mod(5, field.neg_mod(7))));
      dictionary.insert(vec![0, 1, 5], FieldElement::new(&field, field.mul_mod(17, field.neg_mod(7))));
      dictionary.insert(vec![84, 1, 5], FieldElement::new(&field, 40));
      let poly_c = MPolynomial { dictionary };

      assert_eq!(poly_a * poly_b, poly_c);
    }
  }

  #[test]
  fn add () {
    let field = Field::new(FIELD_PRIME);

    {
      let mut dictionary = HashMap::new();
      dictionary.insert(vec![0, 1, 5], FieldElement::new(&field, 17));
      dictionary.insert(vec![5, 23, 0], FieldElement::new(&field, 5));
      let poly_a = MPolynomial { dictionary };

      let mut dictionary = HashMap::new();
      dictionary.insert(vec![42, 0], FieldElement::new(&field, 8));
      dictionary.insert(vec![5, 23], FieldElement::new(&field, 12));
      let poly_b = MPolynomial { dictionary };

      let mut dictionary = HashMap::new();
      dictionary.insert(vec![0, 1, 5], FieldElement::new(&field, 17));
      dictionary.insert(vec![5, 23, 0], FieldElement::new(&field, 17));
      dictionary.insert(vec![42, 0, 0], FieldElement::new(&field, 8));
      let poly_c = MPolynomial { dictionary };

      assert_eq!(poly_a + poly_b, poly_c);
    }
  }

  #[test]
  fn neg () {
    let field = Field::new(FIELD_PRIME);

    {
      let mut dictionary = HashMap::new();
      dictionary.insert(vec![0, 1, 5], FieldElement::new(&field, 17));
      dictionary.insert(vec![5, 23, 0], FieldElement::new(&field, 5));
      let poly_a = MPolynomial { dictionary };

      let mut dictionary = HashMap::new();
      dictionary.insert(vec![0, 1, 5], -FieldElement::new(&field, 17));
      dictionary.insert(vec![5, 23, 0], -FieldElement::new(&field, 5));
      let poly_b = MPolynomial { dictionary };

      assert_eq!(-poly_a, poly_b);
    }
  }

  #[test]
  fn sub () {
    let field = Field::new(FIELD_PRIME);

    {
      let mut dict = HashMap::new();
      dict.insert(vec![0, 1, 5], FieldElement::new(&field, 17));
      dict.insert(vec![5, 23, 0], FieldElement::new(&field, 5));
      let poly_a = MPolynomial::new(dict);

      let mut dict = HashMap::new();
      dict.insert(vec![42, 0], FieldElement::new(&field, 8));
      dict.insert(vec![5, 23], FieldElement::new(&field, 12));
      let poly_b = MPolynomial::new(dict);

      let mut dict = HashMap::new();
      dict.insert(vec![0, 1, 5], FieldElement::new(&field, 17));
      dict.insert(vec![5, 23, 0], FieldElement::new(&field, 5) - FieldElement::new(&field, 12));
      dict.insert(vec![42, 0, 0], -FieldElement::new(&field, 8));
      let poly_c = MPolynomial::new(dict);

      assert_eq!(poly_a - poly_b, poly_c);
    }
  }

  #[test]
  fn variable () {
    let field = Field::new(FIELD_PRIME);

    {
      let mut dict = HashMap::new();
      dict.insert(vec![1, 0, 0], FieldElement::new(&field, 1));
      let a = MPolynomial::new(dict);

      let mut dict = HashMap::new();
      dict.insert(vec![0, 1, 0], FieldElement::new(&field, 1));
      let b = MPolynomial::new(dict);

      let mut dict = HashMap::new();
      dict.insert(vec![0, 0, 1], FieldElement::new(&field, 1));
      let c = MPolynomial::new(dict);

      assert_eq!(MPolynomial::variables(3, &field), vec![a, b, c]);
    }
  }

  #[test]
  fn lift () {
    let field = Field::new(FIELD_PRIME);
    let variables = MPolynomial::variables(4, &field);
    let zero = field.zero();
    let one = field.one();
    let two = FieldElement::new(&field, 2);
    let five = FieldElement::new(&field, 5);

    let upoly = Polynomial::interpolate_domain(&vec![zero, one, two], &vec![two, five, five]);
    let mpoly = MPolynomial::lift(&upoly, 3);

    assert_eq!(upoly.evaluate(&five), mpoly.evaluate(&vec![zero, zero, zero, five]))
  }

  #[test]
  fn evaluate () {
    let field = Field::new(FIELD_PRIME);
    let variables = MPolynomial::variables(4, &field);
    let zero = field.zero();
    let one = field.one();
    let two = FieldElement::new(&field, 2);
    let five = FieldElement::new(&field, 5);

    let mpoly1 = MPolynomial::constant(one)
        * variables[0].clone()
        + MPolynomial::constant(two)
        * variables[1].clone()
        + MPolynomial::constant(five)
        * (variables[2].clone() ^ 3);
    let mpoly2 = MPolynomial::constant(one)
        * variables[0].clone()
        * variables[3].clone()
        + MPolynomial::constant(five)
        * (variables[3].clone() ^ 3)
        + MPolynomial::constant(five);

    let point = vec![zero, five, five, two];

    let eval1 = mpoly1.evaluate(&point);
    let eval2 = mpoly2.evaluate(&point);

    assert_eq!(eval1 * eval2, (mpoly1.clone() * mpoly2.clone()).evaluate(&point));
    assert_eq!(eval1 + eval2, (mpoly1.clone() + mpoly2.clone()).evaluate(&point));
  }

  #[test]
  fn evaluate_symbolic () {
    let field = Field::new(FIELD_PRIME);

    let mpoly = {
      let mut dict = HashMap::new();
      dict.insert(vec![0, 1, 5], FieldElement::new(&field, 17));
      dict.insert(vec![6, 2, 13], FieldElement::new(&field, 8));
      MPolynomial::new(dict)
    };

    let polys = {
      let mut vec = vec![];
      vec.push(Polynomial::new(vec![
        FieldElement::new(&field, 5),
        FieldElement::new(&field, 0),
        FieldElement::new(&field, 2),
      ]));
      vec.push(Polynomial::new(vec![
        FieldElement::new(&field, 2),
        FieldElement::new(&field, 6),
        FieldElement::new(&field, 34),
      ]));
      vec.push(Polynomial::new(vec![
        FieldElement::new(&field, 8),
        FieldElement::new(&field, 9),
        FieldElement::new(&field, 10),
      ]));
      vec
    };
    let poly_res =
      Polynomial::new(vec![FieldElement::new(&field, 17)])
        * (polys[0].clone() ^ 0_u128)
        * (polys[1].clone() ^ 1_u128)
        * (polys[2].clone() ^ 5_u128)
      + Polynomial::new(vec![FieldElement::new(&field, 8)])
        * (polys[0].clone() ^ 6_u128)
        * (polys[1].clone() ^ 2_u128)
        * (polys[2].clone() ^ 13_u128)
    ;

    assert_eq!(mpoly.evaluate_symbolic(&polys), poly_res);
  }

  #[test]
  fn pow () {
    let field = Field::new(FIELD_PRIME);

    let mut dict = HashMap::new();
    dict.insert(vec![1,2,5], FieldElement::new(&field, 3));
    dict.insert(vec![5,3,4], FieldElement::new(&field, 4));
    let mpoly = MPolynomial::new(dict);

    // (3  * x   * y^2 * z^5  + 4  * x^5  * y^3 * z^4) ** 3
    // 1. (3  * x   * y^2 * z^5  + 4  * x^5  * y^3 * z^4)
    // 2. (9  * x^2 * y^4 * z^10 + 24 * x^6  * y^5 * z^9  + 16 * x^10 * y^6 * z^8 )
    //    (27 * x^3 * y^6 * z^15 + 72 * x^7  * y^7 * z^14 + 48 * x^11 * y^8 * z^13) +
    //    (36 * x^7 * y^7 * z^14 + 96 * x^11 * y^8 * z^13 + 64 * x^15 * y^9 * z^12)
    //
    // + 144 * x^11 * y^8 * z^13
    // + 27  * x^3  * y^6 * z^15
    // + 108 * x^7  * y^7 * z^14
    // + 64  * x^15 * y^9 * z^12

    let mut dict = HashMap::new();
    dict.insert(vec![11, 8, 13], FieldElement::new(&field, 144));
    dict.insert(vec![ 3, 6, 15], FieldElement::new(&field,  27));
    dict.insert(vec![ 7, 7, 14], FieldElement::new(&field, 108));
    dict.insert(vec![15, 9, 12], FieldElement::new(&field,  64));

    assert_eq!(mpoly ^ 3, MPolynomial::new(dict));
  }
}

