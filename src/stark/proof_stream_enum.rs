use std::io::Read;
use crate::field::field::Field;
use crate::field::field_element::FieldElement;
use crate::utils::bytes::Bytes;
use crate::utils::digest::{Digest, digest_vec};

#[derive(Debug, Clone, PartialEq)]
pub enum StarkProofStreamEnum<'a> {
    Root(Bytes),
    Codeword(Vec<FieldElement<'a>>),
    Path(Vec<Bytes>),
    Leafs((FieldElement<'a>, FieldElement<'a>, FieldElement<'a>)),
    Value(FieldElement<'a>),
}

impl<'a> StarkProofStreamEnum<'a> {
    pub fn from_bytes(code: u8, b: Bytes, field: &'a Field) -> Self {
        match code {
            0 => StarkProofStreamEnum::Root(b),
            1 => {
                let b = b.bytes();
                StarkProofStreamEnum::Codeword(
                    b
                        .chunks(u128::BITS as usize / 8)
                        .map(|v| u128::from_be_bytes(v.try_into().expect("incorrect size")))
                        .map(|v| FieldElement::new(&field, v))
                        .collect()
                )
            },
            2 => {
                let mut path = vec![];

                let mut b = b.bytes();
                let mut size = [0; 64 / 8];
                while let Ok(_) = b.read_exact(&mut size) {
                    let size = usize::from_be_bytes(size.try_into().expect("incorrect size"));
                    let mut el = vec![0; size];
                    b.read_exact(&mut el).unwrap();

                    path.push(Bytes::new(el));
                }

                StarkProofStreamEnum::Path(path)
            },
            3 => {
                let b = b.bytes();
                let mut iter = b
                    .chunks(u128::BITS as usize / 8)
                    .map(|v| u128::from_be_bytes(v.try_into().expect("incorrect size")))
                    .map(|v| FieldElement::new(&field, v));
                StarkProofStreamEnum::Leafs((
                    iter.next().expect("Leaf to be 3 element"),
                    iter.next().expect("Leaf to be 3 element"),
                    iter.next().expect("Leaf to be 3 element"),
                ))
            },
            4 => {
                let b = b.bytes();
                StarkProofStreamEnum::Value(
                    FieldElement::new(&field, u128::from_be_bytes(b.try_into().expect("incorrect size")))
                )
            },
            _ => panic!("Unknown code"),
        }
    }

    pub fn to_bytes(&self) -> (u8, Bytes, Option<&'a Field>) {
        match self {
            StarkProofStreamEnum::Root(r) => {
                (
                    0,
                    r.clone(),
                    None,
                )
            }
            StarkProofStreamEnum::Codeword(c) => {
                let mut field = None;
                let b = digest_vec(&c, |fe| {
                    if let Some(field) = field {
                        if fe.field != field {
                            panic!("codeword elements in different fields")
                        }
                    } else {
                        field = Some(fe.field)
                    }
                    Bytes::new(fe.value.to_be_bytes().to_vec())
                });
                (
                    1,
                    b,
                    field,
                )
            }
            StarkProofStreamEnum::Path(p) => {
                (
                    2,
                    digest_vec(&p, |b| Bytes::new(b.bytes().len().to_be_bytes().to_vec()) + b.clone()),
                    None,
                )
            }
            StarkProofStreamEnum::Leafs(l) => {
                let field = if l.0.field == l.1.field && l.1.field == l.2.field {
                    Some(l.0.field)
                } else {
                    panic!("leaf elements in different fields")
                };

                let b =
                    Bytes::new(l.0.value.to_be_bytes().to_vec())
                    + Bytes::new(l.1.value.to_be_bytes().to_vec())
                    + Bytes::new(l.2.value.to_be_bytes().to_vec());

                (
                    3,
                    b,
                    field,
                )
            }
            StarkProofStreamEnum::Value(fe) => {
                (
                    4,
                    Bytes::new(fe.value.to_be_bytes().to_vec()),
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

impl Digest for &[StarkProofStreamEnum<'_>] {
    fn digest<'m>(&'m self) -> Bytes {
        let mut field = None;
        let bytes = digest_vec(self, |pse| {
            let (code, bytes, f) = pse.to_bytes();
            if let Some(f) = f {
                if let Some(field) = field {
                    if f != field {
                        panic!("codeword elements in different fields")
                    }
                } else {
                    field = Some(f)
                }
            }

            let bytes = bytes.bytes();
            let bytes_size = bytes.len().to_be_bytes();
            let mut v = Vec::with_capacity(1 + bytes_size.len() + bytes.len());
            v.push(code);
            v.extend(bytes_size);
            v.extend(bytes);

            Bytes::new(v)
        });

        Bytes::new(
            field.map(|fe| fe.order).unwrap_or(0).to_be_bytes().to_vec()
        ) + bytes
    }
}
