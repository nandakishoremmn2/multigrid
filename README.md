Multigrid Solver
=========

Solves laplace equation using multigrid techniques

To run solver - 

`make all && time ./solver_serial`
or
`make cilk && time ./solver_cilk`
or
`make openmp && time ./solver_openmp`

To visualize solution run `python plot.py out.dat`

To plot convergence history run `python plot_data.py data.dat`
