//! Utilities
//!
//! Utility functions and structs for rusty ray tracer
//!
//!
// use tracing::{error, debug}; // , error, info, span, trace, warn};

use rand::Rng;
use rand::distributions::uniform::{SampleRange, SampleUniform};

/// Generate a tuple (x, y) of random values, where x is create within range_x and y in range_y.
///
/// # Example
/// ```rust
/// let (a, b) = random_2d_pos(0..=2, -2..=2);
/// assert!(a <= 2);
/// assert!(a >= 0);
/// assert!(b <= 2);
/// assert!(b >= -2);
/// ```
pub fn random_2d_pos<T, R>(range_x: R, range_y: R) -> (T, T)
where
        T: SampleUniform,
        R: SampleRange<T>
{
    let mut rng = rand::thread_rng();

    let rndm_x = rng.gen_range(range_x);
    let rndm_y = rng.gen_range(range_y);

    (rndm_x, rndm_y)
}

/// Generate a random value between 0 and 1
///
/// # Example
/// ```rust
/// let a = random();
/// assert!(a <= 1);
/// assert!(a >= 0);
/// ```
pub fn random() -> f64 {
    let mut rng = rand::thread_rng();

    rng.gen_range(0.0..1.0)
}

/// Generate a random value between 0 and 2
///
/// # Example
/// ```rust
/// let a = random3();
/// assert!(a <= 2);
/// assert!(a >= 0);
/// ```
pub fn random3() -> u8 {
    let mut rng = rand::thread_rng();

    rng.gen_range(0..=2)
}

#[cfg(test)]
mod tests {
    // use super::*;

    // #[test]
    // fn it_works() {
    //     let result = add(2, 2);
    //     assert_eq!(result, 4);
    // }
}
