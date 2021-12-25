use crate::{ R, C };

#[cfg(feature = "ft_cpu")]
use crate::ft_cpu::MFtStruct;

pub struct Solver {
    shape: Vec<usize>,
    grid_fourier_space: Vec<C>,
    ft_struct: MFtStruct
}
