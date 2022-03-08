# Nonlinear Schrödinger solver (work in progress)

This is a simple numerical solver for the nonlinear Schrödinger equations using a spectral approach. 

# Limitations

This project is in a very early phase, and most features are not implemented yet.

My aim is to design it from scratch as much as possible to maximize flexibility. This will also certainly make it less efficient that solvers built on existing crates like [`num`](https://docs.rs/num/0.4.0/num/) for defining the `Complex` type and [`fftw3`](https://github.com/rust-math/fftw) for the Fourier transforms. If you would like to use this project but performance is an issue, please let me know—I'll be happy to write a version using faster libraries than my custom ones. (Or, of course, feel free to do it yourself and submit a pull request!)

The Fourier transform engine is primarily designed for power-of-two sizes. It should work with non-power-of-two sizes too, but with worse performances.

# Structure

*To write*

# Examples

Running the examples requires a Rust 2021 compiler (see the [rist-lang](https://www.rust-lang.org/) website for installation instructions). Cargo is also recommended. Each example can be run using the command 

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
* test the solver on a 1D quasi-stationary solution
* test the solver on a 2D quasi-stationary solution

## To do

* add documentation on the nonlinear Schrödinger equations
* 1D fft using OpenCL
* multi-dimensional fft using OpenCL
* link the OpenCL fft to the solver
* implement a 1D visualizer (using plotters)
* implement a 2D visualizer (using plotters)
* tests 
* benchmarks
