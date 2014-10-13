#ifndef MULTIGRIDSOLVER_H
#define MULTIGRIDSOLVER_H

#include "multigrid.h"

class MultiGridSolver
{
public:
	/**
		Initialises the solver
		@param grid_density The no. of grid points on the finest grid = ( 2^grid_density + 1 )^2
		@param no_of_grids The no. of grids to use ( < grid_density )
	*/
	MultiGridSolver(int grid_density, int no_of_grids);
	~MultiGridSolver();

	/**
		Solve the equation on the given grid upto a given tolerance level
		@param tol tolerance level
	*/
	void solve(real tol);

private:
	MultiGrid *Grid;

	/**
		V cycle scheme
	*/
	void V(MultiGrid *grid, int v1, int v2);

	/**
		W cycle scheme
	*/
	void W(MultiGrid *grid);

	/**
		FMG scheme
		Depends on V cycle
	*/
	void FMG(MultiGrid *grid);

};

#endif // MULTIGRIDSOLVER_H