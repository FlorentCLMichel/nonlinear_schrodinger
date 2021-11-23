# Nonlinear Schrödinger solver (work in progress)

This is a simple numerical solver for the nonlinear Schrödinger equations.

## To do

* add documentation on the nonlinear Schrödinger equations
* multi-dimensional fft on CPU
* structure for the solver:
    * dimensionality
    * shape of the grid and space step
    * external potential (closure)
    * nonlinear terms (closure)
    * whether to use OpenCL or not
* implement the solver: one time step (take the time step as argument)
* implement the solver: multiple time steps (take the time step and total time as arguments)
* 1D fft using OpenCL
* multi-dimensional fft using OpenCL
* link the OpenCL fft to the solver
* implement a 1D visualizer (using plotters)
* implement a 2D visualizer (using plotters)
* tests 
* benchmarks
