use crate::{ R, C, PI };

#[cfg(feature = "ft_cpu")]
use crate::ft_cpu::MFtStruct;

pub struct Solver {
    shape: Vec<usize>,
    steps: Vec<R>,
    grid_fourier_space: Vec<Vec<R>>,
    potential: Vec<R>,
    ft_struct: MFtStruct
}


impl Solver {

    /// generale a new `solver`
    ///
    /// # Argument 
    ///
    /// * `shape`: shape of the grid
    /// * `steps`: steps in the space directions
    /// * `potential`: array of potential values
    ///
    /// # Error
    ///
    /// Returns a [`SolverError`] if the two arguments have different lengths.
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::{ Solver, R };
    ///
    /// # fn main() -> Result<(), nonlinear_schrodinger::SolverError> {
    /// #
    /// let shape: Vec<usize> = vec![100, 100, 100];
    /// let steps: Vec<R> = vec![0.1, 0.1, 0.1];
    /// let potential: Vec<R> = vec![0.; 1_000_000];
    /// let solver = Solver::new(shape, steps, potential)?;
    /// #
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(shape: Vec<usize>, steps: Vec<R>, potential: Vec<R>) -> Result<Self, SolverError> {

        // check the length of the shape vector
        if shape.len() != steps.len() {
            return Err(SolverError::new(format!(
                    "The 'steps' vector (size {}) must have the same size as the 'shape' one (size {})",
                    steps.len(), shape.len())));
        }

        let ft_struct = MFtStruct::new(shape.clone());
        let grid_fourier_space = generate_fourier_grid(&shape, &steps);

        // check the length of the potential vector
        let n_points = shape.iter().fold(1, |res, a| res * a);
        if potential.len() != n_points {
            return Err(SolverError::new(format!(
                    "The 'potential' vector (size {}) must have the same size as the number of points ({})",
                    potential.len(), n_points)));
        }

        Ok(Solver { shape, steps, grid_fourier_space, potential, ft_struct })
    }
}


#[derive(Debug, Clone)]
pub struct SolverError {
    message: String
}

impl std::fmt::Display for SolverError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", &self.message)
    }
}

impl std::error::Error for SolverError {}

impl SolverError {
    fn new(message: String) -> Self {
        SolverError { message }
    }
}


// assume `shape` and `steps` have the same length
fn generate_fourier_grid(shape: &[usize], steps: &[R]) -> Vec<Vec<R>> {
    let mut grid = Vec::<Vec<R>>::with_capacity(shape.len());
    for (i, &s) in shape.iter().enumerate() {
        let k_max = 2.*PI / steps[i];
        let dk = k_max / (s as R);
        let mut row = Vec::<R>::with_capacity(s);
        for j in 0..=(s/2) {
            row.push((j as R) * dk);
        }
        for j in (s/2+1)..s {
            row.push((j as R) * dk - k_max);
        }
        grid.push(row);
    }
    grid
}
