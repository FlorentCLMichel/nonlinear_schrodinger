# Nonlinear Schrödinger solver (work in progress)

This is a simple numerical solver for the nonlinear Schrödinger equations using a spectral approach. 

# Limitations

My aim is to design it from scratch as much as possible as a learning experience. For this reason, it is less efficient that solvers built on existing well-optimised crates like [`num`](https://docs.rs/num/0.4.0/num/) for defining the `Complex` type and [`fftw3`](https://github.com/rust-math/fftw) for the Fourier transforms. If you would like to use this project but performance is an issue, please let me know—I'll be happy to write a version using faster libraries than my custom ones. (Or, of course, feel free to do it yourself and submit a pull request!)

The Fourier transform engine is primarily designed for power-of-two sizes. It should work with non-power-of-two sizes too, but with worse performances.

# How to use

## Nonlinear Schrödinger solver

The outside-facing part of the crate can be imported using 

```Rust
use nonlinear_schrodinger::prelude::*;
```

This line will import

* `R`: The real type used by the crate (currently `f64`),
* `C`: The complex type used by the crate (currently memory equivalent to `[C; 2]`),
* `Solver`: the solver class,
* `SolverError`: An error type returned by a `Solver` instance when an operation fails.

A new `Solver` instance is built using the `Solver::new` function, with signature

```Rust
pub fn new(shape: Vec<usize>, 
           steps: Vec<R>, 
           potential_fun: Box<dyn Fn(&[R]) -> R>, 
           g: Box<dyn Fn(R) -> R>) 
    -> Result<Solver, SolverError> 
```

It takes four arguments: 

* `shape`: length of the grid in each space direction; for instance, if `shape` is `vec![10, 20, 30]`, the solver will solve the non-linear Schrödinger equation on a 3-dimensional grid with shape 10 by 20 by 30;
* `steps`: space step in each direction (must have the same length as `shape`);
* `potential_fun`: external potential as a function from $\mathbb{R}^N$ to $\mathbb{R}$, where $N$ is the length of `shape`; 
* `g`: nonlinear term, as a function of the density.

A `solver` instance has the following member functions: 

* `pub fn eval_psi(&self, psi_fun: Box<dyn Fn(&[R]) -> C>) -> Vec<C>`: generate a discretized version of the function `psi_fun` on the solver's numerical grid (mostly useful for setting the initial condition), 
* `pub fn mass(&self, psi: &[C]) -> Result<R, SolverError>`: compute the mass from the (discretized) condensate wave function `psi`,
* `pub fn momentum(&mut self, psi: &[C]) -> Result<Vec<R>, SolverError>` : compute the momentum from the (discretized) condensate wave function `psi`,
* `pub fn evolve(&mut self, psi: &mut [C], dt: R, nt: usize) -> Result<(), SolverError>`: evolve the initial configuration `psi` for `nt` time steps of duration `dt` each, 
* `pub fn evolve_i(&mut self, psi: &mut [C], dt: R, nt: usize) -> Result<(), SolverError>`: evolve the initial configuration `psi` for `nt` time steps of duration `dt` each in imaginary time (mostly useful for finding the lowest-energy configuration while preserving the other conserved quantities), 
* `pub fn evolve_oop(&mut self, psi0: &[C], dt: R, nt: usize) -> Result<(), SolverError>`: evolve the initial configuration `psi0` for `nt` time steps of duration `dt` each, without modifying `psi0`, 
* `pub fn evolve_i_oop(&mut self, psi0: &[C], dt: R, nt: usize) -> Result<(), SolverError>`: evolve the initial configuration `psi0` for `nt` time steps of duration `dt` each in imaginary time, without modifying `psi0`.

## Plotters

The crate also provides basic 1D and 2D plot functions, `plotters::plot_1d` and `plotters::plot_2d`, using the `plotters` crate.

### 1D plotter

The function signature is: 

```Rust
pub fn plot_1d(x: &[Vec<R>], 
               y: &[Vec<R>], 
               fname: &str,
               title: &str, 
               xlabel: &str, 
               ylabel: &str, 
               labels: &[String],
               col1: (u8, u8, u8), 
               col2: (u8, u8, u8)) 
    -> Result<(), Box<dyn std::error::Error>> 
```

Arguments: 

* `x`: array of vectors of `x` coordinates (1 for each curve),
* `y`: array of vectors of `y` coordinates (1 for each curve; `y` must have the same length as `x` and, or each `i` between `0` and `x.len()-1`, `y[i]` mush have the same length as `x[i]`),
* `fname`: file name to use for saving the plot, 
* `title`: the plot title,
* `xlabel`: label for the horizontal axis,
* `ylabel`: label for the vertical axis,
* `labels`: label for each curve (must have the same length as `x`), 
* `col1`: colour of the first curve, 
* `col2`: colour of the last curve (if `x.len() > 1`).

### 2D plotter

The function signature is: 

```Rust
pub fn plot_2d(z: &[R], 
               fname: &str, 
               x_min: R, 
               x_max: R, 
               y_min: R, 
               y_max: R,
               nx: usize, 
               ny: usize, 
               show_axes: bool)
    -> Result<(), Box<dyn std::error::Error>> 
```

Arguments: 

* `z`: z coordinate for each point, 
* `fname`: file name to use for saving the plot, 
* `x_min`: minimum value of the x coordinate, 
* `x_max`: maximum value of the x coordinate, 
* `y_min`: minimum value of the y coordinate, 
* `y_max`: maximum value of the y coordinate, 
* `nx`: number of points in the x direction,
* `ny`: number of points in the y direction,
* `show_axes`: whether to show the x and y axes.

# Examples

Running the examples requires a Rust 2021 compiler (see the [rust-lang](https://www.rust-lang.org/) website for installation instructions). Cargo is also recommended. Each example can be run using the command 

```
cargo run --release --example [example]
```

## `1d_harmonic_gaussian`

One-dimensional simulation with a harmonic potential and a Gaussian initial condition. The evolution of the density is saved as an `svg` file (`Results/1d_harmonic_gaussian.svg`).

## `1d_harmonic_imag`

One-dimensional imaginary-time simulation with a harmonic potential and a Gaussian initial condition. The evolution of the density is saved as an `svg` file (`Results/1d_harmonic_imag.svg`).

## `2d_harmonic_gaussian`

Two-dimensional simulation with a symmetric harmonic potential and a Gaussian initial condition. The initial and final densities are saved as `png` files (`Results/2d_harmonic_gaussian_initial_density.png` and `Results/2d_harmonic_gaussian_final_density.png`).

## `2d_harmonic_two_gaussians`

Two-dimensional simulation with a symmetric harmonic potential and an initial condition with two Gaussians. Snapshots of the density profile are saved in the folder `Results/2d_harmonic_two_gaussians`.

# Roadmap

## Done

* custom `Complex` type
* custom 1D Fourier transform
* multi-dimensional fft on CPU
* implementation of a first solver
* implement a 1D visualiser
* implement a 2D visualiser
* test the solver on a 1D quasi-stationary solution
* test the solver on a 2D quasi-stationary solution

## To do

* add documentation on the nonlinear Schrödinger equations
* 1D fft using OpenCL
* multi-dimensional fft using OpenCL
* link the OpenCL fft to the solver
* more tests 
* benchmarks
