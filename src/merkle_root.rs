use serde::Serialize;
use crate::crypto::blake2b512::blake2b512;
use crate::field::field_element::FieldElement;
use crate::utils::bytes::Bytes;

pub struct MerkleRoot;

// codeword - [FieldElement]
// codewords - [[FieldElement]]
impl<'a> MerkleRoot {
  fn commit_ (leafs: &[Bytes]) -> Bytes {
    let len = leafs.len();
    assert_eq!(len & (len - 1), 0, "Leafs len must be power of two");

    if len == 1 {
      leafs.first().unwrap().clone()
    } else {
      let a = MerkleRoot::commit_(&leafs[0..len/2]);
      let b = MerkleRoot::commit_(&leafs[len/2..len]);
      let concat = a + b;
      blake2b512(concat)
    }
  }

  pub fn commit <T>(leafs: &[T]) -> Bytes
    where
      T: Serialize
  {
    let leafs = leafs
      .iter()
      .map(|l| {
        let str = serde_json::to_string(l).unwrap();
        blake2b512(str.as_bytes().into())
      })
      .collect::<Vec<_>>();
    MerkleRoot::commit_(&leafs)
  }

  fn open_ (index: usize, leafs: &[Bytes]) -> Vec<Bytes> {
    let len = leafs.len();
    assert_eq!(len & (len - 1), 0, "length must be power of two");
    assert!(index < len, "cannot open invalid index");

    if len == 2 {
      let str = leafs.get(1 - index).unwrap().clone();
      vec![str]
    } else if index < (len / 2) {
      let mut a = MerkleRoot::open_(index, &leafs[0..len/2]);
      let b = MerkleRoot::commit_(&leafs[len/2..len]);
      a.push(b);
      a
    } else {
      let mut a = MerkleRoot::open_(index - len / 2, &leafs[len/2..len]);
      let b = MerkleRoot::commit_(&leafs[0..len/2]);
      a.push(b);
      a
    }
  }

  pub fn open <T>(index: usize, leafs: &[T]) -> Vec<Bytes>
    where
      T: Serialize
  {
    let leafs = leafs
      .iter()
      .map(|l| {
        let str = serde_json::to_string(l).unwrap();
        blake2b512(str.as_bytes().into())
      })
      .collect::<Vec<_>>();
    MerkleRoot::open_(index, &leafs)
  }

  // leaf = stringified FieldElement<'a>
  fn verify_ (root: &Bytes, index: usize, path: &[Bytes], leaf: Bytes) -> bool {
    let len = path.len();
    assert!(index < (1 << len), "Cannot verify invalid index");

    let path_0 = path.get(0).unwrap().to_owned();
    if len == 1 {
      if index == 0 {
        root == &blake2b512(leaf + path_0)
      } else {
        root == &blake2b512(path_0 + leaf)
      }
    } else {
      if index % 2 == 0 {
        MerkleRoot::verify_(root, index >> 1, &path[1..len], blake2b512(leaf + path_0))
      } else {
        MerkleRoot::verify_(root, index >> 1, &path[1..len], blake2b512(path_0 + leaf))
      }
    }
  }

  pub fn verify <T>(root: &Bytes, index: usize, path: &[Bytes], leaf: &T) -> bool
    where
      T: Serialize
  {
    let str = serde_json::to_string(leaf).unwrap();
    let hash = blake2b512(str.as_bytes().into());
    MerkleRoot::verify_(root, index, path, hash)
  }
}

#[cfg(test)]
mod tests {
  use crate::field::field::{Field, FIELD_PRIME};
  use super::*;

  #[test]
  fn test_commit_one () {
    let field = Field::new(FIELD_PRIME);
    // Serialized to: `{"field":{"prime":270497897142230380135924736767050121217},"value":11}`
    let bytes = MerkleRoot::commit(&[FieldElement {
      field: &field,
      value: 11,
    }]);
    assert_eq!(
      bytes,
      "aeb6f644b81a7b6cde7386b6eefa29952775e6e44ae3368e071884b7129293f3a3ff57d25768525faddd6ace5833b715892afefc470a8a4381b91b6901ae1bf3".into(),
    );
  }

  #[test]
  fn test_commit_two () {
    let field = Field::new(FIELD_PRIME);
    let bytes = MerkleRoot::commit(&[FieldElement {
      field: &field,
      value: 5462,
    }, FieldElement {
      field: &field,
      value: 456,
    }]);
    assert_eq!(
      bytes,
      "6d093f2284637958f1b19981b919fc327d8c35cd30a7411044be22de43bdcd356c76a8965dc80aeff44ad24da5c8fb0a73fd89446d9c2482dd16810600800351".into(),
    );
    let bytes = MerkleRoot::commit(&[FieldElement {
      field: &field,
      value: 652,
    }, FieldElement {
      field: &field,
      value: 23409,
    }]);
    assert_eq!(
      bytes,
      "9bec33e09623d5d30b78140f8094882f64bfc1fbcdd21225ba14d5a0e412b0d6c15e6eb0fe123b2428f5e1750db3bfc60d0910aab8960ede9c7fa5c31ab45f20".into(),
    );
  }

  #[test]
  fn test_commit_four () {
    let field = Field::new(FIELD_PRIME);
    let bytes = MerkleRoot::commit(&[FieldElement {
      field: &field,
      value: 5462,
    }, FieldElement {
      field: &field,
      value: 456,
    }, FieldElement {
      field: &field,
      value: 652,
    }, FieldElement {
      field: &field,
      value: 23409,
    }]);
    assert_eq!(
      bytes,
      "885d536abb872bb162a959541296e86019ae020cd2b3d90ef1bac4259315680bea5052f7cc4bff2dec4c507641748738d5869b9ba2ff9880d7086f430443dfe5".into(),
    );
  }

  #[test]
  fn test_open () {
    let field = Field::new(FIELD_PRIME);
    let bytes = MerkleRoot::open(1, &[FieldElement {
      field: &field,
      value: 5462,
    }, FieldElement {
      field: &field,
      value: 456,
    }, FieldElement {
      field: &field,
      value: 652,
    }, FieldElement {
      field: &field,
      value: 23409,
    }]);
    assert_eq!(bytes, vec![
      "ee9496fd96f3774c94d176fd236055a25a79eff98bf043f12fd61fac7ba41e7dff681ab8fbe6a17fb76cbe4ead02705ee2b749641934e225d909dcd77aed9608".into(),
      "9bec33e09623d5d30b78140f8094882f64bfc1fbcdd21225ba14d5a0e412b0d6c15e6eb0fe123b2428f5e1750db3bfc60d0910aab8960ede9c7fa5c31ab45f20".into(),
    ]);
  }

  #[test]
  fn test_verify () {
    let field = Field::new(FIELD_PRIME);
    let result = MerkleRoot::verify(
      &"885d536abb872bb162a959541296e86019ae020cd2b3d90ef1bac4259315680bea5052f7cc4bff2dec4c507641748738d5869b9ba2ff9880d7086f430443dfe5".into(),
      1,
      &vec![
        "ee9496fd96f3774c94d176fd236055a25a79eff98bf043f12fd61fac7ba41e7dff681ab8fbe6a17fb76cbe4ead02705ee2b749641934e225d909dcd77aed9608".into(),
        "9bec33e09623d5d30b78140f8094882f64bfc1fbcdd21225ba14d5a0e412b0d6c15e6eb0fe123b2428f5e1750db3bfc60d0910aab8960ede9c7fa5c31ab45f20".into(),
      ],
      &FieldElement {
        field: &field,
        value: 456,
      });
    assert!(result, "Root has to be valid");

    let result = MerkleRoot::verify(
      &"885d536abb872bb162a959541296e86019ae020cd2b3d90ef1bac4259315680bea5052f7cc4bff2dec4c507641748738d5869b9ba2ff9880d7086f430443dfe5".into(),
      1,
      &vec![
        "ee9496fd96f3774c94d176fd236055a25a79eff98bf043f12fd61fac7ba41e7dff681ab8fbe6a17fb76cbe4ead02705ee2b749641934e225d909dcd77aed9608".into(),
        "9bec33e09623d5d30b78140f8094882f64bfc1fbcdd21225ba14d5a0e412b0d6c15e6eb0fe123b2428f5e1750db3bfc60d0910aab8960ede9c7fa5c31ab45f20".into(),
      ],
      &FieldElement {
        field: &field,
        value: 5462,
      });
    assert!(!result, "Root has to be invalid because element is invalid");

    let result = MerkleRoot::verify(
      &"885d536abb872bb162a959541296e86019ae020cd2b3d90ef1bac4259315680bea5052f7cc4bff2dec4c507641748738d5869b9ba2ff9880d7086f430443dfe5".into(),
      0,
      &vec![
        "ee9496fd96f3774c94d176fd236055a25a79eff98bf043f12fd61fac7ba41e7dff681ab8fbe6a17fb76cbe4ead02705ee2b749641934e225d909dcd77aed9608".into(),
        "9bec33e09623d5d30b78140f8094882f64bfc1fbcdd21225ba14d5a0e412b0d6c15e6eb0fe123b2428f5e1750db3bfc60d0910aab8960ede9c7fa5c31ab45f20".into(),
      ],
      &FieldElement {
        field: &field,
        value: 456,
      });
    assert!(!result, "Root has to be invalid because index is invalid");
  }
}
