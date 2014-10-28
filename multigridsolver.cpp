// forward declared dependencies
// class MultiGrid;

#include <fstream>
#include <stdio.h>
#include <cstdlib>

#include "multigridsolver.h"

MultiGridSolver::MultiGridSolver(int *data)
{
	int grid_density = data[0];
	int no_of_grids = data[1];

	scheme = data[2];
	switch(scheme)
	{
		case 0:
			v1 		= data[3];
			v2 		= data[4];
			cycles 	= data[5];
			break;
		case 1:
			v1 		= data[3];
			v2 		= data[4];
			mu 		= data[5];
			cycles 	= data[6];
			break;
		case 2:
			v0 		= data[3];
			v1 		= data[4];
			v2 		= data[5];
			break;
	}

	Grid = new MultiGrid(grid_density, no_of_grids-1, 1.0);

	// Boundary conditions applied to the finest grid alone
	Grid->apply_boundary_conditions();
}

MultiGridSolver::~MultiGridSolver()
{
	delete Grid;
}

void MultiGridSolver::solve(real tol)
{
	switch(scheme)
	{
		case 0:
			while(Grid->get_L2norm() > tol)
			{
				V(Grid, v1, v2);
			}
			break;
		case 1:
			while(Grid->get_L2norm() > tol)
			{
				W(Grid, v1, v2, mu);
			}
			break;
		case 2:
			FMG(Grid, v0, v1, v2);
			while(Grid->get_L2norm() > tol)
			{
				V(Grid, v1, v2);
			}
			break;
		default:
			for (int i = 0; i < scheme-2; ++i)
			{
				data.push_back(Grid->relax(1));
			}		
	}
	Grid->save_grid("out.dat");
	save_data("data.dat");
}

void MultiGridSolver::save_data(char *filename)
{
	std::ofstream outfile(filename);

	for (std::vector<GridData>::iterator i = data.begin(); i != data.end(); ++i)
	{
		outfile<<i->iterations<<" ";
		outfile<<i->residue<<" ";
		outfile<<i->wu<<" ";
		outfile<<"\n";
	}

}

void MultiGridSolver::save()
{
	Grid->save_grid("out.dat");
	save_data("data.dat");
}

void MultiGridSolver::V(MultiGrid *grid, int v1, int v2)
{
	MultiGrid *cgrid = grid->grid2; // coarser grid

	data.push_back(grid->relax(v1));
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
	data.push_back(grid->relax(v2));
}

void MultiGridSolver::W(MultiGrid *grid, int v1, int v2, int mu)
{
	MultiGrid *cgrid = grid->grid2; // coarser grid

	data.push_back(grid->relax(v1));
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
	data.push_back(grid->relax(v2));	
}

void MultiGridSolver::FMG(MultiGrid *grid, int v0, int v1, int v2)
{
	MultiGrid *cgrid = grid->grid2; // coarser grid

	 // If coarser grid exists recursively apply FMG()
	if (cgrid)
	{
		// Restrict 'f' to coarse grid
		grid->copy_f_to_temp();
		grid->restrict();
		cgrid->copy_temp_to_f();

		// FMG cycle on coarser grid
		FMG(cgrid, v0, v1, v2);

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

	grid->apply_boundary_conditions();

	// Perform V cycle v0 times more times
	for (int i = 0; i < v0; ++i)
	{
		V(grid, v1, v2);	
	}
}