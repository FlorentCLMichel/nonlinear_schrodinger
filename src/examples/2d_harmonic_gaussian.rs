use nonlinear_schrodinger::prelude::*;
use nonlinear_schrodinger::plotters::plot_2d;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    // parameters
    let shape: Vec<usize> = vec![512, 512];
    let steps: Vec<R> = vec![0.025, 0.025];
    let potential_fun = |x: &[R]| { x[0] * x[0] + x[1] * x[1] };
    let g = |x| x;
    
    // initial condition; should be stationary when the amplitude is small
    let amplitude = 0.5;
    let nx = shape[0];
    let ny = shape[1];
    let dx = steps[0];
    let dy = steps[1];
    let mut psi = Vec::<C>::new();
    for i in 0..nx {
        for j in 0..ny {
            psi.push(C::new(- 2_f64.powf(-0.5) * (
                     ((i as R - (nx/2) as R) * dx).powi(2)
                     + ((j as R - (ny/2) as R) * dy).powi(2)
                   ), 0.).exp()
             * amplitude);
        }
    }
    
    // build the solver
    let mut solver = Solver::new(shape, steps, Box::new(potential_fun), Box::new(g))?;

    // plot the initial density
    plot_2d(
        &psi.iter().map(|e| e.abs2()).collect::<Vec<R>>(), 
        "Results/2d_harmonic_gaussian_initial_density", 
        -(nx as R /2.), nx as R / 2., 
        -(ny as R /2.), ny as R / 2., 
        nx, ny)?;
     
    // evolution
    let dt = 0.001;
    let nt = 100;
    println!("Initial mass: \t{}", solver.mass(&psi)?);
    solver.evolve(&mut psi, dt, nt)?;
    println!("Final mass: \t{}", solver.mass(&psi)?);
    
    // plot the final density
    plot_2d(
        &psi.iter().map(|e| e.abs2()).collect::<Vec<R>>(), 
        "Results/2d_harmonic_gaussian_final_density", 
        -(nx as R /2.), nx as R / 2., 
        -(ny as R /2.), ny as R / 2., 
        nx, ny)?;

    Ok(())
}
