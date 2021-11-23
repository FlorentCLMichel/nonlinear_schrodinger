use crate::C;

/// a structure to perform Fourier transforms
pub struct FtStruct {
    n: usize,
    plan: Vec<usize>,
    butterflies: Vec<usize>, 
    twiddles: Vec<C>,
    twiddles_small_ft: Vec<C>
}

mod plan;
pub use plan::*;
mod fft_error;
pub use fft_error::FFTError;
mod implementation;

mod tests;
