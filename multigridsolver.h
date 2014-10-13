#include "multigrid.h"

class MultiGridSolver
{
public:
	MultiGridSolver(real grid_density, int no_of_grids);
	~MultiGridSolver();

	void solve(real tol);

private:
	MultiGrid *Grid;

	/**
		V cycle scheme
	*/
	void V(MultiGrid *grid);

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