use crate::utils::bytes::Bytes;

/// Replaces serde as a proof of concept to have an optimized size of generated proof.
/// Requires own trait instead of `From/To<String>` because of [Only traits defined in the current crate can be implemented for arbitrary types](https://doc.rust-lang.org/error_codes/E0117.html).
pub trait Digest {
    fn digest<'m>(&'m self) -> Bytes;
}

/// For a plain text serialization:
/// Using two different delimiters results in 5% increase in compressed (gzip) size (from 26.98kb to 28.27kb) (from 27624 characters to 28940)
/// but for now it's negligible and allows a lot easier deserialize implementation.
/// And they also doesn't eliminate the need of brackets (`[`, `(`) because otherwise it's impossible to distinguish array elements.
///
/// Using bytes and hex representation (uncompressed size: 62.08 KB = 63566 characters):
/// Compressed size increased to 29.45kb (30152 characters) in little-endian
/// But in big-endian compressed size decreased to 26.29 KB (26916 characters)
pub fn digest_vec<'m, T, F> (
    vec: &'m [T],
    mut f: F,
) -> Bytes
    where
        F: FnMut(&'m T) -> Bytes,
{
    let mut res = Bytes::new(vec![]);
    let mut iter = vec.iter();
    if let Some(o) = iter.next() {
        res = res + f(o);
        res = iter.fold(res, |res, o| {
            res + f(o)
        })
    }
    res
}
