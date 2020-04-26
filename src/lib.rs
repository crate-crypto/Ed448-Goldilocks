#![forbid(unsafe_code)]
// XXX: Change this to deny later on
#![warn(unused_attributes, unused_imports, unused_mut, unused_must_use)]

// As usual, we will use this file to carefully define the API/ what we expose to the user
mod curve;
pub mod decaf;
pub mod field;
// pub mod ristretto;
pub(crate) mod window;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
