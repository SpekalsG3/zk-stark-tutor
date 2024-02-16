use crate::field::field::Field;
use crate::field::field_element::FieldElement;
use crate::utils::bytes::Bytes;
use crate::utils::stringify::{Stringify, stringify_vec};

#[derive(Debug, Clone, PartialEq)]
pub enum StarkProofStreamEnum<'a> {
    Root(Bytes),
    Codeword(Vec<FieldElement<'a>>),
    Path(Vec<Bytes>),
    Leafs((FieldElement<'a>, FieldElement<'a>, FieldElement<'a>)),
    Value(FieldElement<'a>),
}

impl<'a> StarkProofStreamEnum<'a> {
    pub fn unstringify(str: &str, field: &'a Field) -> Self {
        let mut chars = str.chars();
        match chars.next().expect("nothing to unstringify") {
            '"' => StarkProofStreamEnum::Root(Bytes::unstringify(str)),
            '[' => {
                match chars.next().expect("invalid vector") {
                    ']' => panic!("unknown empty vector"),
                    '"' => {
                        let vec = str[1..str.len()-1]
                            .split(',')
                            .map(|s| Bytes::unstringify(s))
                            .collect();
                        StarkProofStreamEnum::Path(vec)
                    },
                    _ => {
                        let vec = str[1..str.len()-1]
                            .split(',')
                            .map(|s| FieldElement::unstringify(s, &field))
                            .collect();
                        StarkProofStreamEnum::Codeword(vec)
                    },
                }
            },
            '(' => {
                let mut iter = str[1..str.len()-1]
                    .split(',')
                    .map(|s| FieldElement::unstringify(s, &field));
                StarkProofStreamEnum::Leafs((
                    iter.next().expect("has to have 3 elements"),
                    iter.next().expect("has to have 3 elements"),
                    iter.next().expect("has to have 3 elements"),
                ))
            },
            _ => StarkProofStreamEnum::Value(
                FieldElement::unstringify(str, &field)
            )
        }
    }

    pub fn stringify(&self) -> (String, Option<&Field>) {
        match self {
            StarkProofStreamEnum::Root(r) => {
                (
                    format!("{}", r.stringify()),
                    None,
                )
            }
            StarkProofStreamEnum::Codeword(c) => {
                let mut field = None;
                let str = stringify_vec(',', c, |fe| {
                    if let Some(field) = field {
                        if fe.field != field {
                            panic!("codeword elements in different fields")
                        }
                    } else {
                        field = Some(fe.field)
                    }
                    fe.value.to_string()
                });
                (str, field)
            }
            StarkProofStreamEnum::Path(p) => {
                (
                    stringify_vec(',', p, |fe| fe.stringify()),
                    None,
                )
            }
            StarkProofStreamEnum::Leafs(l) => {
                let str = format!(
                    "({},{},{})",
                    l.0.value,
                    l.1.value,
                    l.2.value,
                );
                let field = if l.0.field == l.1.field && l.1.field == l.2.field {
                    Some(l.0.field)
                } else {
                    panic!("leaf elements in different fields")
                };
                (str, field)
            }
            StarkProofStreamEnum::Value(fe) => {
                (
                    format!("{}", fe.value),
                    Some(fe.field),
                )
            }
        }
    }

    pub fn expect_root (self) -> Bytes {
        match self {
            StarkProofStreamEnum::Root(x) => x,
            _ => panic!("expected StarkProofStreamEnum to be Root"),
        }
    }
    pub fn expect_codeword (self) -> Vec<FieldElement<'a>> {
        match self {
            StarkProofStreamEnum::Codeword(x) => x,
            _ => panic!("expected StarkProofStreamEnum to be Codeword"),
        }
    }
    pub fn expect_path (self) -> Vec<Bytes> {
        match self {
            StarkProofStreamEnum::Path(x) => x,
            _ => panic!("expected StarkProofStreamEnum to be Path"),
        }
    }
    pub fn expect_leafs (self) -> (FieldElement<'a>, FieldElement<'a>, FieldElement<'a>) {
        match self {
            StarkProofStreamEnum::Leafs(x) => x,
            _ => panic!("expected StarkProofStreamEnum to be Leafs"),
        }
    }
    pub fn expect_value (self) -> FieldElement<'a> {
        match self {
            StarkProofStreamEnum::Value(x) => x,
            _ => panic!("expected StarkProofStreamEnum to be Value"),
        }
    }
}

impl Stringify for &[StarkProofStreamEnum<'_>] {
    fn stringify<'m>(&'m self) -> String {
        let mut field = None;
        let str = stringify_vec(';', self, |pse| {
            let (str, f) = pse.stringify();
            if let Some(f) = f {
                if let Some(field) = field {
                    if f != field {
                        panic!("codeword elements in different fields")
                    }
                } else {
                    field = Some(f)
                }
            }
            str
        });
        format!("{};{}", field.map(|f| f.to_string()).unwrap_or("_".to_string()), str)
    }
}
