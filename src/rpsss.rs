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
        expansion_factor: usize, // 4
        num_collinearity_checks: usize, // 64
        security_level: usize, // >= 2 * num_collinearity_checks
        transition_constraints_degree: usize, // 2
    ) -> Self {
        let rp = RescuePrime::new(field, 2, 1, security_level, 27);
        Self {
            stark: Stark::new(
                field,
                expansion_factor,
                num_collinearity_checks,
                security_level,
                rp.m,
                rp.N + 1,
                transition_constraints_degree,
            ),
            rp,
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
        sps: SignatureProofStream<'m>,
    ) -> Result<(), String> {
        let boundary_constraints = self.rp.boundary_constraints(output_element);
        let transition_constraints = self.rp.transition_constraints(self.stark.omicron);
        self.stark.verify::<SignatureProofStream>(&transition_constraints, &boundary_constraints, sps)
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

    pub fn verify<'m, T: Into<Bytes>>(&'m self, pk: FieldElement, document: T, signature: Bytes ) -> Result<(), String> {
        let mut sps = SignatureProofStream::new(document);
        sps.ps = self.stark.deser_independent_proof_stream(signature);
        self.stark_verify(pk, sps)
    }
}

#[cfg(test)]
mod tests {
    use std::time::SystemTime;
    use crate::field::field::{Field, FIELD_PRIME};
    use crate::rpsss::RPSSS;
    use crate::utils::bytes::Bytes;

    #[test]
    fn rpsss () {
        let field = Field::new(FIELD_PRIME);
        let rpsss = RPSSS::new(&field, 4, 64, 128, 2);

        let time = SystemTime::now();
        println!("Started at {}", time.duration_since(std::time::UNIX_EPOCH).unwrap().as_millis());

        let time = SystemTime::now();
        let (sk, pk) = rpsss.keygen();
        let time = SystemTime::now().duration_since(time).unwrap();
        println!("Generated keys in {}ms", time.as_millis());

        let doc = Bytes::from("Hello, World!".as_bytes());
        let time = SystemTime::now();
        let signature = rpsss.sign(sk, doc.clone());
        let time = SystemTime::now().duration_since(time).unwrap();
        assert!(signature.is_ok(), "Failed to generate signature - {}", signature.unwrap_err());
        println!("Generated signature in {}ms", time.as_millis());

        let signature = signature.unwrap();
        let time = SystemTime::now();
        let is_verified = rpsss.verify(pk, doc, signature.clone());
        let time = SystemTime::now().duration_since(time).unwrap();
        assert!(is_verified.is_ok(), "Failed to verify signature for valid doc - {}", is_verified.unwrap_err());
        println!("Verified signature for valid doc in {}ms", time.as_millis());

        let doc = Bytes::from("Malicious document".as_bytes());
        let time = SystemTime::now();
        let is_verified = rpsss.verify(pk, doc, signature.clone());
        let time = SystemTime::now().duration_since(time).unwrap();
        assert!(is_verified.is_err(), "Shouldn't successfully verified signature for invalid doc");
        println!("Verified signature for invalid doc in {}ms", time.as_millis());

        println!("\tsize of signature in bytes: {}\n\tin kb: {}", signature.bytes().len(), signature.bytes().len() / (2_usize.pow(13)));
    }
}
