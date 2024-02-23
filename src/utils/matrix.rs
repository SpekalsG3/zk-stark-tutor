use crate::field::field::Field;
use crate::field::field_element::FieldElement;

// source - https://rosettacode.org/wiki/Reduced_row_echelon_form
pub fn rref(matrix: &mut Vec<Vec<FieldElement<'_>>>) {
    let mut lead = 0;
    let row_count = matrix.len();
    let column_count = matrix[0].len();
    'outer: for r in 0..row_count {
        if column_count <= lead {
            break
        }

        let mut i = r;
        while matrix[i][lead].is_zero() {
            i += 1;
            if row_count == i {
                i = r;
                lead += 1;
                if column_count == lead {
                    break 'outer;
                }
            }
        }

        { // Swap rows i and r
            let tmp = matrix[i].clone();
            matrix[i] = matrix[r].clone();
            matrix[r] = tmp;
        }

        if !matrix[r][lead].is_zero() {
            matrix[r] = matrix[r]
                .clone()
                .into_iter()
                .map(|el| el / matrix[r][lead])
                .collect();
        }
        for i in 0..row_count {
            if i != r {
                let hold = matrix[i][lead];
                for k in 0..column_count {
                    // Subtract M[i, lead] multiplied by row r from row i
                    matrix[i][k] = matrix[i][k] - hold * matrix[r][k];
                }
            }
        }
        lead += 1;
    }
}

// source - https://rosettacode.org/wiki/Matrix_transposition
pub fn transpose<'m>(matrix: Vec<Vec<FieldElement<'m>>>) -> Vec<Vec<FieldElement<'m>>> {
    let m = matrix.len();
    let n = matrix[0].len();

    (0..n)
        .map(|col| {
            (0..m)
                .map(|row| {
                    matrix[row][col]
                })
                .collect()
        })
        .collect()
}

pub fn inverse<'m>(
    field: &'m Field,
    matrix: Vec<Vec<FieldElement<'m>>>,
) -> Result<Vec<Vec<FieldElement<'m>>>, String> {
    let n = matrix.len();
    let mut matrix = matrix
        .into_iter()
        .enumerate()
        .map(|(i, row)| {
            if row.len() == n {
                let mut v = vec![field.zero(); n];
                v[i] = field.one();
                Ok(row.into_iter().chain(v).collect())
            } else {
                Err("Inverse exists only for square matrices")
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    rref(&mut matrix);

    let (identity, matrix) = matrix
        .into_iter()
        .fold((Vec::with_capacity(n), Vec::with_capacity(n)), |mut acc, row| {
            acc.0.push(row[..n].to_vec());
            acc.1.push(row[n..].to_vec());
            acc
        });

    let one = field.one();
    for (i, row) in identity.iter().enumerate() {
        if
            row[..i].iter().position(|el| !el.is_zero()).is_some()
            || &row[i] != &one
            || row[(i+1)..].iter().position(|el| !el.is_zero()).is_some()
        {
            return Err("Couldnt construct identity matrix to find inverse".to_string())
        }
    }

    Ok(matrix)
}

#[cfg(test)]
mod tests {
    use crate::field::field::{Field, FIELD_PRIME};
    use crate::field::field_element::FieldElement;
    use crate::utils::matrix::{inverse, rref, transpose};

    #[test]
    fn test_transpose () {
        let field = Field::new(FIELD_PRIME);
        let matrix = vec![
            vec![FieldElement::new(&field, 1), FieldElement::new(&field,  2), FieldElement::new(&field,  3), FieldElement::new(&field,  4)],
            vec![FieldElement::new(&field, 5), FieldElement::new(&field,  6), FieldElement::new(&field,  7), FieldElement::new(&field,  8)],
            vec![FieldElement::new(&field, 9), FieldElement::new(&field, 10), FieldElement::new(&field, 11), FieldElement::new(&field, 12)],
        ];
        let res = transpose(matrix);
        let check = vec![
            vec![FieldElement::new(&field, 1), FieldElement::new(&field, 5), FieldElement::new(&field,  9)],
            vec![FieldElement::new(&field, 2), FieldElement::new(&field, 6), FieldElement::new(&field, 10)],
            vec![FieldElement::new(&field, 3), FieldElement::new(&field, 7), FieldElement::new(&field, 11)],
            vec![FieldElement::new(&field, 4), FieldElement::new(&field, 8), FieldElement::new(&field, 12)],
        ];
        assert_eq!(res, check);
    }

    #[test]
    fn test_rref () {
        let field = Field::new(FIELD_PRIME);
        let mut matrix = vec![
            vec![ FieldElement::new(&field, 1), FieldElement::new(&field, 2), -FieldElement::new(&field, 1), -FieldElement::new(&field, 4)],
            vec![ FieldElement::new(&field, 2), FieldElement::new(&field, 3), -FieldElement::new(&field, 1), -FieldElement::new(&field, 11)],
            vec![-FieldElement::new(&field, 2), FieldElement::new(&field, 0), -FieldElement::new(&field, 3),  FieldElement::new(&field, 22)],
        ];
        rref(&mut matrix);

        let res = vec![
            vec![FieldElement::new(&field, 1), FieldElement::new(&field, 0), FieldElement::new(&field, 0), -FieldElement::new(&field, 8)],
            vec![FieldElement::new(&field, 0), FieldElement::new(&field, 1), FieldElement::new(&field, 0),  FieldElement::new(&field, 1)],
            vec![FieldElement::new(&field, 0), FieldElement::new(&field, 0), FieldElement::new(&field, 1), -FieldElement::new(&field, 2)],
        ];
        assert_eq!(matrix, res);
    }

    #[test]
    fn test_inverse () {
        let field = Field::new(FIELD_PRIME);

        {
            let matrix = vec![
                vec![FieldElement::new(&field, 1), FieldElement::new(&field, 2), FieldElement::new(&field, 3)],
                vec![FieldElement::new(&field, 0), FieldElement::new(&field, 1), FieldElement::new(&field, 4)],
                vec![FieldElement::new(&field, 5), FieldElement::new(&field, 6), FieldElement::new(&field, 0)],
            ];
            let matrix = inverse(&field, matrix);
            assert!(matrix.is_ok());
            assert_eq!(matrix.unwrap(), vec![
                vec![-FieldElement::new(&field, 24), FieldElement::new(&field, 18), FieldElement::new(&field, 5)],
                vec![FieldElement::new(&field, 20), -FieldElement::new(&field, 15), -FieldElement::new(&field, 4)],
                vec![-FieldElement::new(&field, 5), FieldElement::new(&field, 4), FieldElement::new(&field, 1)],
            ]);
        }
        {
            let matrix = vec![
                vec![FieldElement::new(&field, 270497897142230380135924736767050121214), FieldElement::new(&field,  4)],
                vec![FieldElement::new(&field, 270497897142230380135924736767050121205), FieldElement::new(&field, 13)],
            ];
            let matrix = inverse(&field, matrix);
            assert!(matrix.is_ok());
            assert_eq!(matrix.unwrap(), vec![
                vec![FieldElement::new(&field, 210387253332845851216830350818816760948), FieldElement::new(&field, 60110643809384528919094385948233360270)],
                vec![FieldElement::new(&field, 90165965714076793378641578922350040407), FieldElement::new(&field, 180331931428153586757283157844700080811)],
            ])
        }
    }
}
