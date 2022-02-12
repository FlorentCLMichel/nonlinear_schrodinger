# Nonlinear Schrödinger solver (work in progress)

This is a simple numerical solver for the nonlinear Schrödinger equations using a spectral approach. 

# Limitations

This project is in a very early phase, and most features are not implemented yet.

My aim is to design it from scratch as much as possible to maximize flexibility. This will also certainly make it less efficient that solvers built on existing crates like [`num`](https://docs.rs/num/0.4.0/num/) for defining the `Complex` type and [`fftw3`](https://github.com/rust-math/fftw) for the Fourier transforms. If you would like to use this project but performance is an issue, please let me know—I'll be happy to write a version using faster libraries than my custom ones. (Or, of course, feel free to do it yourself and submit a pull request!)

The Fourier transform engine is primarily designed for power-of-two sizes. It should work with non-power-of-two sizes too, but with significanty worse performances, especially for sizes which can not be factorized into products of small integers.

# Roadmap

## Done

* custom `Complex` type
* custom 1D Fourier transform
* multi-dimensional fft on CPU
* implementation of a first solver

## To do

* test the solver on 1D gray solitons
* test the solver on 2D configurations
* add documentation on the nonlinear Schrödinger equations
* 1D fft using OpenCL
* multi-dimensional fft using OpenCL
* link the OpenCL fft to the solver
* implement a 1D visualizer (using plotters)
* implement a 2D visualizer (using plotters)
* tests 
* benchmarks
