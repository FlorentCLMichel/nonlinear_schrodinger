use crate::C;

/// a structure to perform one-dimensional Fourier transforms
pub struct FtStruct {
    n: usize,
    plan: Vec<usize>,
    butterflies: Vec<usize>, 
    twiddles: Vec<C>,
    twiddles_small_ft: Vec<C>
}

/// a structure to perform multi-dimensional Fourier transforms
pub struct MFtStruct {
    dimension: usize, 
    total_length: usize, 
    shape: Vec<usize>,
    ft_structs: Vec<FtStruct>
}

mod plan;
pub use plan::*;
mod fft_error;
pub use fft_error::FFTError;
mod implementation;

mod tests;
