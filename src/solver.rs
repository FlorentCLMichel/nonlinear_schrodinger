use crate::{ R, C, PI };

#[cfg(feature = "ft_cpu")]
use crate::ft_cpu::{ MFtStruct, FFTError };

pub struct Solver {
    shape: Vec<usize>,
    steps: Vec<R>,
    n_points: usize,
    dimensions: Vec<R>,
    infinitesimal_volume: R,
    infinitesimal_volume_k: R,
    grid: Vec<Vec<R>>,
    grid_k: Vec<Vec<R>>,
    kinetic: Vec<R>,
    potential: Vec<R>,
    ft_struct: MFtStruct
}


impl Solver {

    /// generate a new `solver`
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
        let n_points = shape.iter().product();

        // check the length of the shape vector
        if shape.len() != steps.len() {
            return Err(SolverError::new(format!(
                    "The 'steps' vector (size {}) must have the same size as the 'shape' one (size {})",
                    steps.len(), shape.len())));
        }
        
        let dimensions: Vec<R> = steps.iter().zip(shape.iter()).map(|(x, y)| x * (*y as R)).collect();
        let infinitesimal_volume = steps.iter().fold(1., |res, a| res * a);
        let infinitesimal_volume_k = infinitesimal_volume / n_points as R;

        let ft_struct = MFtStruct::new(shape.clone());

        // generate the kinetic and potential terms
        let mut grid = Vec::<Vec<R>>::with_capacity(n_points);
        let mut grid_k = Vec::<Vec<R>>::with_capacity(n_points);
        let mut kinetic = Vec::<R>::with_capacity(n_points);
        let mut potential = Vec::<R>::with_capacity(n_points);
        let grid_fourier_space = generate_fourier_grid(&shape, &steps);
        let mut coordinates: Vec<R> = vec![0.; shape.len()];
        let mut coordinates_k: Vec<R> = vec![0.; shape.len()];
        for i in 0..n_points {
            let mut k = i;
            let mut k2: R = 0.;
            for j in 0..shape.len() {
                let l = k % shape[j];
                let kj = grid_fourier_space[j][l];
                coordinates_k[j] = kj;
                k2 += kj * kj;
                coordinates[j] = 0.5 * steps[j] * ((2*l + 1) as R - shape[j] as R);
                k /= shape[j];
            }
            grid.push(coordinates.clone());
            potential.push(potential_fun(&coordinates));
            grid_k.push(coordinates_k.clone());
            kinetic.push(0.5*k2);
        }

        Ok(Solver { shape, steps, n_points, dimensions, infinitesimal_volume, 
                    infinitesimal_volume_k, grid, grid_k, kinetic, potential, ft_struct })
    }
    

    /// compute the array representing a wave function from a closure
    ///
    /// # Argument
    ///
    /// * `psi`: A closure of taking an array of `dimension` numbers of type `R` and returning a
    /// number of type `C`.
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::{ Solver, R, C };
    ///
    /// # fn main() -> Result<(), nonlinear_schrodinger::SolverError> {
    /// #
    ///
    /// // build the solver structure
    /// let shape: Vec<usize> = vec![100, 100, 100];
    /// let steps: Vec<R> = vec![0.1, 0.1, 0.1];
    /// let potential_fun = |x: &[R]| { x[0] + x[1] + x[2] };
    /// let solver = Solver::new(shape, steps, Box::new(potential_fun))?;
    ///
    /// // example of wave function
    /// let psi_fun = |x: &[R]| { C::new(x[0].cos(), x[0].sin()) };
    /// let psi = solver.eval_psi(Box::new(psi_fun)); 
    /// #
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval_psi(&self, psi_fun: Box<dyn Fn(&[R]) -> C>) -> Vec<C> {
        self.grid.iter().map(|x| psi_fun(x)).collect()
    }


    /// compute the mass
    ///
    /// # Argument
    ///
    /// * `psi`: A complex array with length `self.n_points` representing the wave function.
    ///
    /// # Return 
    ///
    /// If `psi` has the right number of points, return `Ok(mass)`. 
    /// 
    /// Otherwise, return a `SolverError`.
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::{ Solver, R, C };
    ///
    /// # fn main() -> Result<(), nonlinear_schrodinger::SolverError> {
    /// #
    ///
    /// // build the solver structure
    /// let shape: Vec<usize> = vec![100, 100, 100];
    /// let steps: Vec<R> = vec![0.1, 0.1, 0.1];
    /// let potential_fun = |x: &[R]| { x[0] + x[1] + x[2] };
    /// let solver = Solver::new(shape, steps, Box::new(potential_fun))?;
    ///
    /// // example of wave function
    /// let psi_fun = |x: &[R]| C::new(5.*x[0].cos(), 5.*x[0].sin());
    /// let psi = solver.eval_psi(Box::new(psi_fun));
    ///
    /// // compute the mass
    /// let mass = solver.mass(&psi)?;
    /// assert!((mass - 25_000.).abs() < 1e-11);
    /// #
    /// # Ok(())
    /// # }
    /// ```
    pub fn mass(&self, psi: &[C]) -> Result<R, SolverError> {

        if psi.len() != self.n_points {
            return Err(SolverError { message: 
                format!("Wrong number of points for the input array in the computation of the mass: expected {}, got {}",
                        self.n_points, psi.len())
            });
        }

        let mut result: R = 0.;
        for x in psi {
            result += x.abs2();
        }
        result *= self.infinitesimal_volume;

        Ok(result)
    }
    

    /// compute the momentum
    ///
    /// # Argument
    ///
    /// * `psi`: A complex array with length `self.n_points` representing the wave function.
    ///
    /// # Return 
    ///
    /// If `psi` has the right number of points, return `Ok(momentum)`, where `momentum`. 
    /// 
    /// Otherwise, return a `SolverError`.
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::{ Solver, R, C, PI };
    ///
    /// # fn main() -> Result<(), nonlinear_schrodinger::SolverError> {
    /// #
    ///
    /// // build the solver structure
    /// let shape: Vec<usize> = vec![100, 100, 100];
    /// let steps: Vec<R> = vec![0.1, 0.1, 0.1];
    /// let potential_fun = |x: &[R]| { x[0] + x[1] + x[2] };
    /// let solver = Solver::new(shape, steps, Box::new(potential_fun))?;
    ///
    /// // example of wave function
    /// let psi_fun = |x: &[R]| C::new(5.*(0.2*x[0]*PI).cos(), 5.*(0.2*x[0]*PI).sin());
    /// let psi = solver.eval_psi(Box::new(psi_fun));
    ///
    /// // compute the momentum
    /// let momentum = solver.momentum(&psi)?;
    /// assert!((momentum[0] - 5_000.*PI).abs() < 1e-11);
    /// assert!((momentum[1] - 0.).abs() < 1e-11);
    /// assert!((momentum[2] - 0.).abs() < 1e-11);
    /// #
    /// # Ok(())
    /// # }
    /// ```
    pub fn momentum(&self, psi: &[C]) -> Result<Vec<R>, SolverError> {

        if psi.len() != self.n_points {
            return Err(SolverError { message: 
                format!("Wrong number of points for the input array in the computation of the momentum: expected {}, got {}",
                        self.n_points, psi.len())
            });
        }

        let mut result: Vec<R> = vec![0.; self.shape.len()];
        let psi_k = self.ft_struct.ft(&psi)?;
        for (k, psi) in self.grid_k.iter().zip(psi_k.iter()) {
            let psi2 = psi.abs2();
            for i in 0..self.shape.len() {
                result[i] += k[i] * psi2 * self.infinitesimal_volume_k;
            }
        }

        Ok(result)
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

impl std::convert::From<FFTError> for SolverError {
    fn from(e: FFTError) -> Self {
        SolverError { message: format!("FFTError: {}", e) }
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
