#![forbid(unsafe_code)]
#![deny(unused_attributes, unused_imports, unused_mut, unused_must_use)]

// As usual, we will use this file to carefully define the API/ what we expose to the user

mod field;
#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
