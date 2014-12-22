Multigrid Solver
=========

Solves laplace equation using multigrid techniques

To run solver - 

`make all && time ./solver_serial`

or

in branch cilk -
`make cilk && time ./solver_cilk`

or

in branch openmp -
`make openmp && time ./solver_openmp`

or

in branch cuda -
`make cuda && time ./solver_cuda`

To visualize solution run `python plot.py out.dat`

To plot convergence history run `python plot_data.py data.dat`
