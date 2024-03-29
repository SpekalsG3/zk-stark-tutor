use crate::crypto::blake2b512::blake2b512;
use crate::crypto::shake256::shake256;
use crate::proof_stream::{IndependentProofStream, ProofStream};
use crate::stark::proof_stream_enum::StarkProofStreamEnum;
use crate::utils::bytes::Bytes;
use crate::utils::digest::Digest;

#[derive(Debug)]
pub struct SignatureProofStream<'a> {
    pub(crate) ps: IndependentProofStream<StarkProofStreamEnum<'a>>,
    pub(crate) prefix: Bytes,
}

impl<'a> SignatureProofStream<'a> {
    pub fn new<T: Into<Bytes>>(
        document: T,
    ) -> Self {
        Self {
            ps: IndependentProofStream::new(),
            prefix: blake2b512(document.into()),
        }
    }

    fn digest_prefix(&self) -> Bytes {
        let size = self.prefix.bytes().len().to_be_bytes();
        Bytes::new(size.to_vec()) + self.prefix.clone()
    }
}

impl<'a> ProofStream<StarkProofStreamEnum<'a>> for SignatureProofStream<'a> {
    fn digest(&self) -> Bytes {
        self.ps.digest()
    }

    fn fiat_shamir_prover(
        &self,
        num_bytes: usize,
    ) -> Bytes {
        let bytes = self.digest();

        shake256(self.digest_prefix() + bytes, num_bytes)
    }

    fn fiat_shamir_verifier(
        &self,
        num_bytes: usize,
    ) -> Bytes {
        let slice = &self.ps.objects[0..self.ps.read_index];
        let bytes = slice.digest();

        shake256(self.digest_prefix() + bytes, num_bytes)
    }

    fn push (&mut self, obj: StarkProofStreamEnum<'a>) {
        self.ps.push(obj)
    }

    fn pull (&mut self) -> Option<StarkProofStreamEnum<'a>> {
        self.ps.pull()
    }

}
