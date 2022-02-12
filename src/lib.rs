/// generic complex number type
pub mod complex;
use crate::complex::Complex;

/// type of real nuber to use
pub type R = f64;

/// specific type of complex number to use
pub type C = Complex<R>;

/// constant pi (ratio of the circumpherence of a circle to its diameter in Euclidean geometry)
pub const PI: R = std::f64::consts::PI;

#[cfg(feature = "ft_cpu")]
mod ft_cpu;
#[cfg(feature = "ft_cpu")]
pub use ft_cpu::*;

mod solver;
pub use solver::*;

pub mod plotters;
pub mod prelude;
