#ifndef MULTIGRIDSOLVER_H
#define MULTIGRIDSOLVER_H

#include <vector>
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

	/**
		Save the data variable to file
		@param filename data file name
	*/
	void save_data(char *filename);

private:
	MultiGrid *Grid;

	/**
		Holds the residue, iteration and work units spend data
	*/
	std::vector<GridData> data;

	/**
		V cycle scheme
	*/
	void V(MultiGrid *grid, int v1, int v2);

	/**
		W cycle scheme
	*/
	void W(MultiGrid *grid, int v1, int v2, int mu);

	/**
		FMG scheme
		Depends on V cycle
	*/
	void FMG(MultiGrid *grid, int v0, int v1, int v2);

};

#endif // MULTIGRIDSOLVER_H