use std::fmt::Write;

/// Replaces serde as a proof of concept to have an optimized size of generated proof.
/// Requires own trait instead of `From/To<String>` because of [Only traits defined in the current crate can be implemented for arbitrary types](https://doc.rust-lang.org/error_codes/E0117.html).
pub trait Stringify {
    fn stringify<'m>(&'m self) -> String;
}

/// Using two different delimiters results in 5% increase in compressed size (from 26.98kb to 28.27kb) (from 27624 characters to 28940)
/// but for now it's negligible and allows a lot easier deserialize implementation.
/// And they also doesn't eliminate the need of brackets (`[`, `(`) because otherwise it's impossible to distinguish array elements.
pub fn stringify_vec<'m, T, F: FnMut(&'m T) -> String> (
    delimiter: char,
    vec: &'m [T],
    mut f: F,
) -> String {
    let mut str = String::from("[");
    let mut iter = vec.iter();
    if let Some(o) = iter.next() {
        str.write_fmt(format_args!("{}", f(o))).unwrap();
        iter.for_each(|o| {
            str.write_fmt(format_args!("{}{}", delimiter, f(o))).unwrap();
        })
    }
    str + "]"
}
