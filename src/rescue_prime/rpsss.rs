use std::io::Read;
use rand::{RngCore, thread_rng};
use crate::field::field::Field;
use crate::field::field_element::FieldElement;
use crate::rescue_prime::proof_stream::SignatureProofStream;
use crate::rescue_prime::rescue_prime::RescuePrime;
use crate::stark::stark::Stark;
use crate::utils::bytes::Bytes;

// Rescue-Prime Stark Signature Scheme
pub struct RPSSS<'a> {
    rp: RescuePrime<'a>,
    stark: Stark<'a>,
}

impl<'a> RPSSS<'a> {
    pub fn new (
        field: &'a Field,
        expansion_factor: usize,
        num_collinearity_checks: usize,
        security_level: usize,
        transition_constraints_degree: usize,
    ) -> Self {
        let rp = RescuePrime::new(field, 2, 1, 27);
        Self {
            stark: Stark::new(
                field,
                expansion_factor, // 4
                num_collinearity_checks, // 64
                security_level, // 2 * num_colinearity_checks
                rp.m,
                rp.N + 1,
                transition_constraints_degree, // 3
            ),
            rp,
        }
    }

    pub fn deser_signature_proof_stream<'m> (&'m self, proof: Bytes) -> SignatureProofStream<'m> {
        let mut b = proof.bytes();

        let prefix = {
            let mut prefix_size = [0; 64 / 8];
            b.read_exact(&mut prefix_size).unwrap();
            let prefix_size = usize::from_be_bytes(prefix_size);

            let mut prefix = vec![0; prefix_size];
            b.read_exact(&mut prefix).unwrap();
            Bytes::new(prefix)
        };

        SignatureProofStream {
            ps: self.stark.deser_default_proof_stream(proof),
            prefix,
        }
    }

    fn stark_prove<'m> (
        &'m self,
        input_elements: FieldElement<'m>,
        proof_stream: SignatureProofStream<'m>,
    ) -> Result<Bytes, String> {
        let output_element = self.rp.hash(input_elements.clone());
        let trace = self.rp.trace(input_elements);

        let transition_constraints = self.rp.transition_constraints(self.stark.omicron);
        let boundary_constraints = self.rp.boundary_constraints(output_element);
        self.stark.prove(trace, &transition_constraints, &boundary_constraints, proof_stream)
    }

    fn stark_verify<'m>(
        &'m self,
        output_element: FieldElement<'m>,
        stark_proof: Bytes,
    ) -> Result<(), String> {
        let boundary_constraints = self.rp.boundary_constraints(output_element);
        let transition_constraints = self.rp.transition_constraints(self.stark.omicron);
        let proof_stream = self.deser_signature_proof_stream(stark_proof);
        self.stark.verify::<SignatureProofStream>(&transition_constraints, &boundary_constraints, proof_stream)
    }

    pub fn keygen<'m>(&'m self) -> (FieldElement<'m>, FieldElement<'m>) {
        let mut thread_rng = thread_rng();
        let mut bytes = vec![0; 17]; // todo: shouldn't 17 be self.num_colinearity_tests?
        thread_rng.fill_bytes(&mut bytes);
        let sk = self.stark.field.sample(&Bytes::from(bytes));
        let pk = self.rp.hash(sk);
        (sk, pk)
    }

    pub fn sign<'m, T: Into<Bytes>>(&'m self, sk: FieldElement<'m>, document: T) -> Result<Bytes, String> {
        let sps = SignatureProofStream::new(document);
        self.stark_prove(sk, sps)
    }

    pub fn verify(&self, pk: Vec<FieldElement>, signature: Bytes ) -> Result<(), String> {
        assert_eq!(pk.len(), 1, "Currently don't support more than 1 input element");
        self.stark_verify(pk[0], signature)
    }
}

