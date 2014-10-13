#include <fstream>
#include "multigrid.h"
#include <cmath>
#include <cstdlib>

/*
	Constructors and destructors
*/
MultiGrid::MultiGrid(int grid_density, int no_of_child_grids)
{
	n = pow(2, grid_density) + 1;

	h = pow(.5, grid_density);
	h2 = h*h;

	// Allocate memory for the grid
	v = initialise();
	f = initialise();
	temp = initialise();

	// Initialise the coarser child grids
	grid2 = ( no_of_child_grids > 0 ) ? new MultiGrid(grid_density-1, no_of_child_grids-1) : NULL;
}

MultiGrid::~MultiGrid()
{
	// Free memory of variables
	deallocate(v);
	deallocate(f);
	deallocate(temp);

	// Destroy child grids
	delete grid2;
}

/*
	Private methods
*/
real **MultiGrid::initialise()
{
	real **var = new real*[n];
	for (int i = 0; i < n; ++i)
	{
		var[i] = new real[n];
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			var[i][j] = 0;
		}
	}
	return var;
}

void MultiGrid::deallocate(real **var)
{
	for (int i = 0; i < n; ++i)
	{
		delete [] var[i];
	}
	delete [] var;
}

void MultiGrid::copy(real **src, real **dest)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			dest[i][j] = src[i][j];
		}
	}
}

/*
	Public methods
*/
void MultiGrid::relax(int v)
{
	for (int k = 0; k < v; ++k)
	{
		/* Gauss, Jacobi, red-black, etc.... */
	}
}

int MultiGrid::getSize()
{
	return n;
}

void MultiGrid::interp()
{
	// Interpolate values from *grid2.temp to temp
}

void MultiGrid::restrict()
{
	// Restrict values from temp to *grid2.temp
}

void MultiGrid::apply_boundary_conditions()
{
	for (int i = 0; i < n; ++i)
	{
		v[i][0] = (real)i*h2/n;
		v[0][i] = (real)i*h2/n;
		v[i][n-1] = 1.0 - (real)i*h2/n;
		v[n-1][0] = 1.0 - (real)i*h2/n;
	}
}

void MultiGrid::set_v(real val)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			v[i][j] = val;
		}
	}
}

void MultiGrid::calc_res_to_temp()
{
	for (int i = 1; i < n-1; ++i)
	{
		for (int j = 1; j < n-1; ++j)
		{
			// r = f - Av
			temp[i][j] = f[i][j] - ( 4*v[i][j] - ( v[i+1][j+1] + v[i-1][j-1] + v[i-1][j+1] + v[i+1][j-1] ) )/h2;
		}
	}

	// Set residual value at boundary to zero
	for (int i = 0; i < n; ++i)
	{
		temp[0][i] = temp[i][0] = temp[n-1][i] = temp[i][n-1] = 0;
	}
}

void MultiGrid::copy_temp_to_f()
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			f[i][j] = temp[i][j];
		}
	}
}

void MultiGrid::copy_v_to_temp()
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			temp[i][j] = v[i][j];
		}
	}
}

void MultiGrid::add_temp_to_v()
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			v[i][j] = v[i][j] + temp[i][j];
		}
	}
}

void MultiGrid::save_grid(char *filename)
{
	std::ofstream outfile(filename);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			outfile<<v[i][j]<<" ";
		}
		outfile<<"\n";
	}
}