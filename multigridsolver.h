#ifndef MULTIGRIDSOLVER_H
#define MULTIGRIDSOLVER_H

#include <vector>
#include "multigrid.h"

class MultiGridSolver
{
public:
	/**
		Initialises the solver
		@param data contains the scheme to use and the parameters
	*/
	MultiGridSolver(int *data);
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

	/**
		Save the grid and convergence data and exits
		@param signum
	*/
	void save();

private:
	MultiGrid *Grid;

	/**
		The scheme to use (0, V cycle) (1, Mu cycle) (2, FMG)
	*/
	int scheme;

	/**
		Scheme parameters
	*/
	int v0, v1, v2, mu, cycles;

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