use crate::{ R, C, PI };

#[cfg(feature = "ft_cpu")]
use crate::ft_cpu::MFtStruct;

pub struct Solver {
    shape: Vec<usize>,
    steps: Vec<R>,
    kinetic: Vec<R>,
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
    /// let potential_fun = |x: &[R]| { x[0] + x[1] + x[2] };
    /// let solver = Solver::new(shape, steps, Box::new(potential_fun))?;
    /// #
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(shape: Vec<usize>, steps: Vec<R>, potential_fun: Box<dyn Fn(&[R]) -> R>) 
        -> Result<Self, SolverError> 
    {
        let n_points = shape.iter().fold(1, |res, a| res * a);

        // check the length of the shape vector
        if shape.len() != steps.len() {
            return Err(SolverError::new(format!(
                    "The 'steps' vector (size {}) must have the same size as the 'shape' one (size {})",
                    steps.len(), shape.len())));
        }

        let ft_struct = MFtStruct::new(shape.clone());

        // generate the kinetic and potential terms
        let mut kinetic = Vec::<R>::with_capacity(n_points);
        let mut potential = Vec::<R>::with_capacity(n_points);
        let grid_fourier_space = generate_fourier_grid(&shape, &steps);
        let mut coordinates: Vec<R> = vec![0.; shape.len()];
        for i in 0..n_points {
            let mut k = i;
            let mut k2: R = 0.;
            for j in 0..shape.len() {
                let l = k % shape[j];
                let kj = grid_fourier_space[j][l];
                k2 += kj * kj;
                coordinates[j] = 0.5 * steps[j] * ((2*l + 1) as R - shape[j] as R);
                k /= shape[j];
            }
            potential.push(potential_fun(&coordinates));
            kinetic.push(0.5*k2);
        }

        Ok(Solver { shape, steps, kinetic, potential, ft_struct })
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


#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn new_solver_1() {

        // create the solver
        let shape: Vec<usize> = vec![301, 201, 101];
        let steps: Vec<R> = vec![0.1, 0.1, 0.1];
        let potential_fun = |x: &[R]| { x[0] + 2.*x[1] + 3.*x[2] };
        let solver = Solver::new(shape, steps, Box::new(potential_fun)).unwrap();

        assert_eq!(solver.kinetic[0], 0.);
        assert!((solver.kinetic[1] - 0.021786965709185015).abs() < 1e-16);
        assert!((solver.kinetic[2] - 0.08714786283674006).abs() < 1e-16);
        assert!((solver.kinetic[301] - 0.04885821836632439).abs() < 1e-16);
        assert!((solver.kinetic[301*201] - 0.1935026840719411).abs() < 1e-16);
        assert!((solver.kinetic[1*301*201+1*301+1] - 0.2641478681474505).abs() < 1e-16);
        assert_eq!(solver.potential[150+100*301+50*201*301], 0.);
        assert!((solver.potential[151+100*301+50*201*301] - 0.1).abs() < 1e-16);
        assert!((solver.potential[150+101*301+50*201*301] - 0.2).abs() < 1e-16);
        assert!((solver.potential[150+100*301+51*201*301] - 0.3).abs() < 1e-16);
        assert!((solver.potential[149+100*301+50*201*301] + 0.1).abs() < 1e-16);
        assert!((solver.potential[150+99*301+50*201*301] + 0.2).abs() < 1e-16);
        assert!((solver.potential[150+100*301+49*201*301] + 0.3).abs() < 1e-16);
    }
}
