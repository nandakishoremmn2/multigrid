#include <fstream>
#include "multigrid.h"
#include <cmath>
#include <cstdlib>

/*
	Constructors and destructors
*/
MultiGrid::MultiGrid(int grid_density, int no_of_child_grids, float wu)
{
	n = pow(2, grid_density) + 1;

	WU = wu;

	h = pow(.5, grid_density);
	h2 = h*h;

	// Allocate memory for the grid
	v = initialise();
	f = initialise();
	temp = initialise();

	// Initialise the coarser child grids
	grid2 = ( no_of_child_grids > 0 ) ? new MultiGrid(grid_density-1, no_of_child_grids-1, WU/4.0) : NULL;
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

real MultiGrid::norm2(real **data)
{
	// L-squared norm
	real norm2val = 0;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			norm2val += data[i][j]*data[i][j]*h2;
		}
	}

	return sqrt(norm2val);
}

/*
	Public methods
*/
GridData MultiGrid::relax(int vn)
{
	for (int k = 0; k < vn; ++k)
	{
		for (int i = 1; i < n-1; ++i)
		{
			for (int j = 1; j < n-1; ++j)
			{
				v[i][j] = ( v[i][j+1] + v[i][j-1] + v[i+1][j] + v[i-1][j] + h2*f[i][j])/4.0;
			}
		}
	}

	calc_res_to_temp();

	GridData data = {
		norm2(temp),
		vn*WU,
		vn
	};

	return data;
}

int MultiGrid::getSize()
{
	return n;
}

void MultiGrid::interp()
{
	// Interpolate values from *grid2.temp to temp
	real **temp2 = grid2->temp;
	int n2 = grid2->getSize();

	for (int i = 0; i < n2-1; ++i)
	{
		for (int j = 0; j < n2-1; ++j)
		{
			// Between 4 points
			temp[2*i+1][2*j+1] = ( temp2[i][j] + temp2[i][j+1] + temp2[i+1][j] + temp2[i+1][j+1] ) / 4.0;
		}
	}

	for (int i = 0; i < n2-1; ++i)
	{
		for (int j = 0; j < n2-1; ++j)
		{
			// Between 2 points
			temp[2*i+1][2*j] = ( temp2[i][j] + temp2[i+1][j] ) / 2.0;
			temp[2*i][2*j+1] = ( temp2[i][j] + temp2[i][j+1] ) / 2.0;
		}
	}

	for (int i = 1; i < n2-1; ++i)
	{
		for (int j = 1; j < n2-1; ++j)
		{
			temp[2*i][2*j] = temp2[i][j];
		}
	}
}

void MultiGrid::restrict()
{
	// Restrict values from temp to *grid2.temp
	real **temp2 = grid2->temp;
	int n2 = grid2->getSize();

	for (int i = 1; i < n2-1; ++i)
	{
		for (int j = 1; j < n2-1; ++j)
		{
			// Weighted average from all 9 points
			temp2[i][j] = 4.0*temp[2*i][2*j];
			temp2[i][j] += 2.0*( temp[2*i+1][2*j] + temp[2*i-1][2*j] + temp[2*i][2*j-1] + temp[2*i][2*j+1] );
			temp2[i][j] += 1.0*( temp[2*i+1][2*j+1] + temp[2*i-1][2*j-1] + temp[2*i+1][2*j-1] + temp[2*i-1][2*j+1] );
			temp2[i][j] /= 16.0;
		}
	}

}

void MultiGrid::apply_boundary_conditions()
{
	for (int i = 0; i < n; ++i)
	{
		v[0][i] = v[i][0] = (real)i*i*h2;
		// v[0][i] = (real)i*i*h2;
		v[n-1][i] = v[i][n-1] = 1.0 - (real)i*i*h2;
		// v[n-1][i] = 1.0 - (real)i*i*h2;
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
			temp[i][j] = f[i][j] - ( 4*v[i][j] - ( v[i+1][j] + v[i-1][j] + v[i][j+1] + v[i][j-1] ) )/h2;
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
	copy(temp, f);
}

void MultiGrid::copy_v_to_temp()
{
	copy(v, temp);
}

void MultiGrid::add_temp_to_v()
{
	for (int i = 1; i < n-1; ++i)
	{
		for (int j = 1; j < n-1; ++j)
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