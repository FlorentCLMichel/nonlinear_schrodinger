//! finding the ground state using imaginary-time evolution

use nonlinear_schrodinger::prelude::*;
use nonlinear_schrodinger::plotters::plot_1d;

const COL1: (u8, u8, u8) = (0, 0, 255);
const COL2: (u8, u8, u8) = (0, 255, 0);

fn main() -> Result<(), Box<dyn std::error::Error>> {

    // parameters
    let shape: Vec<usize> = vec![512];
    let steps: Vec<R> = vec![0.025];
    let potential_fun = |x: &[R]| { x[0] * x[0] };
    let g = |x| x;
    
    // initial condition
    let amplitude = 0.1;
    let sigma = 2.;
    let nx = shape[0];
    let dx = steps[0];
    let mut psi = (0..nx).map(|i| 
                              C::new(- 2_f64.powf(-0.5) * ((i as R - (nx/2) as R) * dx / sigma).powi(2), 0.).exp()
                              * amplitude / sigma.powf(0.5)
                          )
                         .collect::<Vec<C>>();
    
    // build the solver
    let mut solver = Solver::new(shape, steps, Box::new(potential_fun), Box::new(g))?;
     
    // evolution
    let n_steps = 10;
    let dt = 0.002;
    let nt = 100;
    let mut densities = vec![
        psi.iter().map(|x| x.abs2()).collect::<Vec<R>>()
    ];
    for _ in 0..n_steps {
        solver.evolve_i(&mut psi, dt, nt)?;
        densities.push(psi.iter().map(|x| x.abs2()).collect::<Vec<R>>());
    }
    
    // plot the densities
    plot_1d(&vec![solver.grid.iter().map(|p| p[0]).collect()], 
            &densities, 
            "Results/1d_harmonic_imag",
            "Density", "x", "ρ", 
            &(0..=n_steps).map(|i| format!("t = {:.2}i", ((i*nt) as R) * dt)).collect::<Vec<String>>(),
            COL1, COL2)?;

    Ok(())
}
