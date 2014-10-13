// forward declared dependencies
// class MultiGrid;

#include "multigridsolver.h"

MultiGridSolver::MultiGridSolver(int grid_density, int no_of_grids)
{
	Grid = new MultiGrid(grid_density, no_of_grids-1);
}

MultiGridSolver::~MultiGridSolver()
{
	delete Grid;
}

void MultiGridSolver::solve(real tol)
{
	V(Grid);
}

void MultiGridSolver::V(MultiGrid *grid)
{
	int v1 = 1;
	int v2 = 3;

	grid->relax(v1);
	if (grid->grid2)
	{
		grid->calc_res_to_temp();
		grid->restrict();
		grid->grid2->copy_temp_to_f();
		grid->grid2->set_v(0);
		V(grid->grid2);
		grid->grid2->copy_v_to_temp();
		grid->interp();
		grid->add_temp_to_v();
	}
	grid->relax(v2);
}

void MultiGridSolver::W(MultiGrid *grid)
{
	
}

void MultiGridSolver::FMG(MultiGrid *grid)
{
	
}