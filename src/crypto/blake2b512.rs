use blake2::{Blake2b512, Digest};
use crate::utils::bytes::Bytes;

pub fn blake2b512 (bytes: Bytes) -> Bytes {
  let mut hasher = Blake2b512::new();

  hasher.update(bytes.bytes());

  Bytes::new(
    hasher
    .finalize()
    .to_vec()
  )
}

#[cfg(test)]
mod tests {
  use crate::crypto::blake2b512::blake2b512;
  use crate::utils::bytes::Bytes;

  #[test]
  fn test () {
    assert_eq!(
      blake2b512(Bytes::from(vec![0u8])).to_hex(),
      "2fa3f686df876995167e7c2e5d74c4c7b6e48f8068fe0e44208344d480f7904c36963e44115fe3eb2a3ac8694c28bcb4f5a0f3276f2e79487d8219057a506e4b"
    );
    assert_eq!(
      blake2b512(Bytes::from(vec![0u8, 0])).to_hex(),
      "5ba7f7e4ade7e5803c59d184326420823f7f860effcfba0bb896d568f59b8d85181cfff25929d40b18e01069c2ef5c31754f1d821a1f3f80f896f4dde374a2f1"
    );
  }
}
