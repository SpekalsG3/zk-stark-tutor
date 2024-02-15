use std::io::Read;
use sha3::{digest::{Update, ExtendableOutput}, Shake256};
use crate::utils::bytes::Bytes;

pub const PROOF_BYTES: usize = 32;

pub fn shake256 (bytes: Bytes, num_bytes: usize) -> Bytes {
    let mut hasher = Shake256::default();

    hasher.update(bytes.bytes());

    let mut buf = vec![0u8; num_bytes];
    let reader = hasher.finalize_xof();
    reader
        .take(num_bytes as u64)
        .read(&mut buf)
        .unwrap(); // under the hood it asserts so Result is useless
    Bytes::new(buf)
}
