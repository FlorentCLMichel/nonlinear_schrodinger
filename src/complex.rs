//! A simple implementation of complex numbers
//!
//! Provide a `Complex` struct made of two real numbers and basic operations.

use std::ops;

mod real_number;
use real_number::RealNumber;

/// The complex type
///
/// `T` must be an ‘acceptable’ real number type
#[derive(Copy, Clone)]
pub struct Complex<T> 
where T: RealNumber
{
    pub real: T,
    pub imag: T
}

mod complex_impl;


#[allow(dead_code)]
/// root mean-square-error between two complex arrays
///
/// # Arguments
///
/// * `a`: first input
/// * `b`: second input
///
/// If the two arrays have different lengths, the smaller one is implicitely padded with 0s.
///
/// # Example
///
/// ```
/// use nonlinear_schrodinger::complex::{ Complex, rmse };
///
/// # fn main() {
/// #
/// let a = vec![Complex::new(3.,0.), Complex::new(0.,2.)];
/// let b = vec![Complex::new(0.,0.), Complex::new(0.,-2.)];
///
/// let dist = rmse(&a, &b);
///
/// assert_eq!(dist, 5.);
/// #
/// # }
/// ```
pub fn rmse<'a, T: RealNumber>(mut a: &'a [Complex<T>], mut b: &'a [Complex<T>]) -> T {
    let mut res = T::from(0);

    // exchange the arrays if a is larger than b
    if a.len() > b.len() {
        std::mem::swap(&mut a, &mut b);
    }

    // compute the rmse
    for i in 0..a.len() {
        res += (a[i] - b[i]).abs2();
    }
    for e in b.iter().skip(a.len()) {
        res += e.abs2();
    }
    res.sqrt()
}


#[cfg(test)]
mod tests;
