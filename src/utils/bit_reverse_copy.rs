use crate::utils::bit_iter::BitIter;

pub fn bit_reverse_copy<T: Clone>(zero: T, inputs: Vec<T>) -> Vec<T> {
    let orig_n = inputs.len();
    let n = orig_n.next_power_of_two();
    let mut new_inputs = vec![zero.clone(); n];

    // pad the input with zeros if the input is not the size of the power two
    let inputs = inputs
        .into_iter()
        .chain(
            vec![zero.clone(); n - orig_n]
        );

    // bit-reverse-copy
    let biggest_bit = (n - 1).ilog2() as usize;
    for (k, el) in inputs.enumerate() {
        let k_rev = {
            let mut n = 0;
            let bit_iter = BitIter(Some(biggest_bit), k);
            for (i, bit) in bit_iter.enumerate() {
                n = n + (if bit { 1 } else { 0 } << i)
            }
            n
        };
        new_inputs[k_rev] = el;
    }

    new_inputs
}
