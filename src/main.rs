use crate::field::field::Field;

mod crypto {
    pub mod blake2b512;
}
mod field {
    pub mod field;
    pub mod field_element;
    pub mod polynomial;
}
mod utils {
    pub mod gcd;
    pub mod xgcd;
    pub mod bit_iter;
    pub mod bytes;
}
pub mod fri;
pub mod m_polynomial;
pub mod merkle_root;
pub mod proof_stream;

fn main() {
    let p = 1 + 407 * (1 << 119);
    let field = Field::new(p);
}
