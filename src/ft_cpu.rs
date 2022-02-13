use crate::{ R, C };
use std::cell::RefCell;

/// a structure to perform one-dimensional Fourier transforms
pub struct FtStruct {
    n: usize,
    butterflies: Vec<usize>, 
    twiddles: Vec<C>,
    c_bar: Vec<C>,
    f_c_tilde: Vec<C>,
    buffer1: RefCell<Vec<C>>,
    buffer2: RefCell<Vec<C>>,
    is_power_2: bool,
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
