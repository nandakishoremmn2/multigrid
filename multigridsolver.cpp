// forward declared dependencies
// class MultiGrid;

#include "multigridsolver.h"

MultiGridSolver::MultiGridSolver(int grid_density, int no_of_grids)
{
	Grid = new MultiGrid(grid_density, no_of_grids-1);

	// Boundary conditions applied to the finest grid alone
	Grid->apply_boundary_conditions();
}

MultiGridSolver::~MultiGridSolver()
{
	delete Grid;
}

void MultiGridSolver::solve(real tol)
{
	V(Grid, 1, 10);
	// Grid->save_grid("out.dat");
}

void MultiGridSolver::V(MultiGrid *grid, int v1, int v2)
{
	MultiGrid *cgrid = grid->grid2; // coarser grid

	grid->relax(v1);
	 // If coarser grid exists recursively apply V()
	if (cgrid)
	{
		// Restrict residue to coarse grid
		grid->calc_res_to_temp();
		grid->restrict();
		cgrid->copy_temp_to_f();

		// Set error 0 on coarse grid
		cgrid->set_v(0);

		// Continue V cycle
		V(cgrid, v1, v2);

		// Interpolate and add error to fine grid
		cgrid->copy_v_to_temp();
		grid->interp();
		grid->add_temp_to_v();
	}
	// Relax fine grid few more times
	grid->relax(v2);
}

void MultiGridSolver::W(MultiGrid *grid)
{
	
}

void MultiGridSolver::FMG(MultiGrid *grid)
{
	
}