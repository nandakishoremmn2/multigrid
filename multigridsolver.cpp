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
	// V(Grid, 1, 5);
	// W(Grid, 1, 5, 2);
	FMG(Grid, 1);
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

void MultiGridSolver::W(MultiGrid *grid, int v1, int v2, int mu)
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

		// Continue M cycle 'mu' times
		for (int i = 0; i < mu; ++i)
		{
			W(cgrid, v1, v2, mu);
		}

		// Interpolate and add error to fine grid
		cgrid->copy_v_to_temp();
		grid->interp();
		grid->add_temp_to_v();
	}
	// Relax fine grid few more times
	grid->relax(v2);	
}

void MultiGridSolver::FMG(MultiGrid *grid, int v0)
{
	MultiGrid *cgrid = grid->grid2; // coarser grid

	 // If coarser grid exists recursively apply V()
	if (cgrid)
	{
		// Restrict residue to coarse grid
		grid->calc_res_to_temp();
		grid->restrict();
		cgrid->copy_temp_to_f();

		// FMG cycle on coarser grid
		FMG(cgrid, v0);

		// Interpolate solution ('v') to fine grid
		cgrid->copy_v_to_temp();
		grid->interp();
		grid->set_v(0);
		grid->add_temp_to_v();
	}
	else // Base case
	{
		grid->set_v(0);
	}

	// Perform V cycle v0 times more times
	for (int i = 0; i < v0; ++i)
	{
		V(grid, 1, 5);	
	}
}