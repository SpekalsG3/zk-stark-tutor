use std::fmt::Write;

pub trait Stringify {
    fn stringify<'m>(&'m self) -> String;
}

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
