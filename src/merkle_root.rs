use crate::crypto::blake2b512::blake2b512;
use crate::utils::bytes::Bytes;

pub struct MerkleRoot;

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
      T: Into<Bytes> + Clone
  {
    let leafs = leafs
      .iter()
      .map(|l| {
        blake2b512(l.clone().into())
      })
      .collect::<Vec<_>>();
    MerkleRoot::commit_(&leafs)
  }

  fn open_ (index: usize, leafs: &[Bytes]) -> Vec<Bytes> {
    let len = leafs.len();
    assert_eq!(len & (len - 1), 0, "length must be power of two");
    assert!(index < len, "cannot open invalid index");

    if len == 2 {
      let commit = leafs.get(1 - index).unwrap().clone();
      vec![commit]
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
      T: Into<Bytes> + Clone
  {
    let leafs = leafs
      .iter()
      .map(|l| {
        blake2b512(l.clone().into())
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
        T: Into<Bytes> + Clone
  {
    let hash = blake2b512(leaf.clone().into());
    MerkleRoot::verify_(root, index, path, hash)
  }
}

#[cfg(test)]
mod tests {
  use crate::field::{
    field_element::FieldElement,
    field::{Field, FIELD_PRIME},
  };
  use super::*;

  #[test]
  fn commit_one () {
    let field = Field::new(FIELD_PRIME);
    // Serialized to: `{"field":{"prime":270497897142230380135924736767050121217},"value":11}`

    let bytes = MerkleRoot::commit(&[FieldElement {
      field: &field,
      value: 11,
    }]);
    assert_eq!(
      bytes,
      "7aa7e388f8145d395ac616bb526eaa35b10069f49e2b36d7327157d1d4af360dfbbfea805aa7e405ed025ce5eadd56c27c40b92991727a5a16b51df5604ad006".into(),
    );

    let bytes = MerkleRoot::commit(&[FieldElement {
      field: &field,
      value: 5462,
    }]);
    assert_eq!(
      bytes.to_hex(),
      "1f069c52b4f26c7714dbd9babacbff542d1333190e3246dec47ee9f30bb649046406f3e0ae8f4cafd52bc1a1305061b451a8746ad3ad240c2524a82a3fcd28c0"
    );
  }

  #[test]
  fn commit_two () {
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
      "e79bb3f920912c56d27de11b3aaedf523d75877d7ec34d7b5819142ba69ce421e665b176fbbbd7b81e90dce61b1f629830eec87c3f7d0644c412af12f47548fe".into(),
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
      "9b70e42c4b3aea3efddaeda6c1883b38c8969e40ca17566d612156c0457961e7c30d811e2adefd941da7b5329d24ecf015dcffb3e39e379dc988564d588a2341".into(),
    );
  }

  #[test]
  fn commit_four () {
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
      "b36f5edab7ea2100fc298d9811bf1a745745282e80243e3a919e71ef6c30f690606b445557ad7843d3251c8e92b83b584d94b738334ffa7d88babd6e47471ac5".into(),
    );
  }

  #[test]
  fn open () {
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
      "1f069c52b4f26c7714dbd9babacbff542d1333190e3246dec47ee9f30bb649046406f3e0ae8f4cafd52bc1a1305061b451a8746ad3ad240c2524a82a3fcd28c0".into(),
      "9b70e42c4b3aea3efddaeda6c1883b38c8969e40ca17566d612156c0457961e7c30d811e2adefd941da7b5329d24ecf015dcffb3e39e379dc988564d588a2341".into(),
    ]);
  }

  #[test]
  fn verify () {
    let field = Field::new(FIELD_PRIME);
    let result = MerkleRoot::verify(
      &"b36f5edab7ea2100fc298d9811bf1a745745282e80243e3a919e71ef6c30f690606b445557ad7843d3251c8e92b83b584d94b738334ffa7d88babd6e47471ac5".into(),
      1,
      &vec![
        "1f069c52b4f26c7714dbd9babacbff542d1333190e3246dec47ee9f30bb649046406f3e0ae8f4cafd52bc1a1305061b451a8746ad3ad240c2524a82a3fcd28c0".into(),
        "9b70e42c4b3aea3efddaeda6c1883b38c8969e40ca17566d612156c0457961e7c30d811e2adefd941da7b5329d24ecf015dcffb3e39e379dc988564d588a2341".into(),
      ],
      &FieldElement {
        field: &field,
        value: 456,
      });
    assert!(result, "Root has to be valid");

    let result = MerkleRoot::verify(
      &"b36f5edab7ea2100fc298d9811bf1a745745282e80243e3a919e71ef6c30f690606b445557ad7843d3251c8e92b83b584d94b738334ffa7d88babd6e47471ac5".into(),
      1,
      &vec![
        "1f069c52b4f26c7714dbd9babacbff542d1333190e3246dec47ee9f30bb649046406f3e0ae8f4cafd52bc1a1305061b451a8746ad3ad240c2524a82a3fcd28c0".into(),
        "9b70e42c4b3aea3efddaeda6c1883b38c8969e40ca17566d612156c0457961e7c30d811e2adefd941da7b5329d24ecf015dcffb3e39e379dc988564d588a2341".into(),
      ],
      &FieldElement {
        field: &field,
        value: 5462,
      });
    assert!(!result, "Root has to be invalid because element is invalid");

    let result = MerkleRoot::verify(
      &"b36f5edab7ea2100fc298d9811bf1a745745282e80243e3a919e71ef6c30f690606b445557ad7843d3251c8e92b83b584d94b738334ffa7d88babd6e47471ac5".into(),
      0,
      &vec![
        "1f069c52b4f26c7714dbd9babacbff542d1333190e3246dec47ee9f30bb649046406f3e0ae8f4cafd52bc1a1305061b451a8746ad3ad240c2524a82a3fcd28c0".into(),
        "9b70e42c4b3aea3efddaeda6c1883b38c8969e40ca17566d612156c0457961e7c30d811e2adefd941da7b5329d24ecf015dcffb3e39e379dc988564d588a2341".into(),
      ],
      &FieldElement {
        field: &field,
        value: 456,
      });
    assert!(!result, "Root has to be invalid because index is invalid");
  }
}
