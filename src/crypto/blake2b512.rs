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
