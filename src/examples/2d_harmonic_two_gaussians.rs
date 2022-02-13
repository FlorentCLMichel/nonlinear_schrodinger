use nonlinear_schrodinger::prelude::*;
use nonlinear_schrodinger::plotters::plot_2d;

const FOLDER_NAME: &str = "Results/2d_harmonic_two_gaussians";

fn main() -> Result<(), Box<dyn std::error::Error>> {

    // create the folder to save the results in
    std::fs::create_dir_all(FOLDER_NAME)?;

    // parameters
    let shape: Vec<usize> = vec![512, 512];
    let steps: Vec<R> = vec![0.025, 0.025];
    let potential_fun = |x: &[R]| { x[0] * x[0] + x[1] * x[1] };
    let g = |x| x;
    
    // initial condition
    let amplitude = 0.1;
    let position = 3.;
    let nx = shape[0];
    let ny = shape[1];
    let dx = steps[0];
    let dy = steps[1];
    let mut psi = Vec::<C>::new();
    for i in 0..nx {
        for j in 0..ny {
            psi.push(
                C::new(- 2_f64.powf(-0.5) * (
                     ((i as R - (nx/2) as R) * dx - position).powi(2)
                     + ((j as R - (ny/2) as R) * dy).powi(2)
                ), 0.).exp() * amplitude
                + 
                C::new(- 2_f64.powf(-0.5) * (
                     ((i as R - (nx/2) as R) * dx + position).powi(2)
                     + ((j as R - (ny/2) as R) * dy).powi(2)
                ), 0.).exp() * amplitude);
        }
    }
    
    // build the solver
    let mut solver = Solver::new(shape, steps, Box::new(potential_fun), Box::new(g))?;

    // plot the initial density
    plot_2d(
        &psi.iter().map(|e| e.abs2()).collect::<Vec<R>>(), 
        &format!("{}/0.0", FOLDER_NAME), 
        -(nx as R /2.)*dx, (nx as R / 2.)*dx, 
        -(ny as R /2.)*dx, (ny as R / 2.)*dy, 
        nx, ny, false)?;
     
    // evolution
    let n_steps = 20;
    let dt = 0.001;
    let nt = 100;
    println!("Mass at time 0.0: \t{}", solver.mass(&psi)?);
    for i in 0..n_steps {
        solver.evolve(&mut psi, dt, nt)?;
        println!("Mass at time {:.1}: \t{}", (((i+1)*nt) as R)*dt, solver.mass(&psi)?);
        
        // plot the density
        plot_2d(
            &psi.iter().map(|e| e.abs2()).collect::<Vec<R>>(), 
            &format!("{}/{:.1}", FOLDER_NAME, (((i+1)*nt) as R)*dt), 
            -(nx as R /2.)*dx, (nx as R / 2.)*dx, 
            -(ny as R /2.)*dx, (ny as R / 2.)*dy, 
            nx, ny, false)?;
    }

    Ok(())
}
