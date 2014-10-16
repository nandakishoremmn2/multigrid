#include <fstream>
#include "multigrid.h"
#include <cmath>
#include <cstdlib>
#include <cilk/cilk.h>

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

	var[0:n][0:n] = 0;

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
	dest[0:n][0:n] = src[0:n][0:n];
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
	real **vt = v, **ft = f;
	for (int k = 0; k < vn; ++k)
	{
		vt[1:n-2][1:n-2] = ( vt[1:n-2][2:n-2] + vt[1:n-2][0:n-2] + vt[2:n-2][1:n-2] + vt[0:n-2][1:n-2] + h2*ft[1:n-2][1:n-2])/4.0;
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
	real **temp2 = grid2->temp, **tmp = temp;
	int n2 = grid2->getSize();

	// Between 4 points
	tmp[1:n2-1:2][1:n2-1:2] = ( temp2[0:n2-1][0:n2-1] + temp2[0:n2-1][1:n2-1] + temp2[1:n2-1][0:n2-1] + temp2[1:n2-1][1:n2-1] ) / 4.0;

	// Between 2 points
	tmp[1:n2-1:2][2:n2-2:2] = ( temp2[0:n2-1][1:n2-2] + temp2[1:n2-1][1:n2-2] ) / 2.0;
	tmp[2:n2-2:2][1:n2-1:2] = ( temp2[1:n2-2][0:n2-1] + temp2[1:n2-2][1:n2-1] ) / 2.0;

	tmp[2:n2-2:2][2:n2-2:2] = temp2[1:n2-2][1:n2-2];
}

void MultiGrid::restrict()
{
	// Restrict values from temp to *grid2.temp
	real **temp2 = grid2->temp, **tmp = temp;
	int n2 = grid2->getSize();

	// Weighted average from all 9 points
	temp2[1:n2-2][1:n2-2] = 4.0*tmp[2:n2-2:2][2:n2-2:2];
	temp2[1:n2-2][1:n2-2] += 2.0*( tmp[3:n2-2:2][2:n2-2:2] + tmp[1:n2-2:2][2:n2-2:2] + tmp[2:n2-2:2][1:n2-2:2] + tmp[2:n2-2:2][3:n2-2:2] );
	temp2[1:n2-2][1:n2-2] += 1.0*( tmp[3:n2-2:2][3:n2-2:2] + tmp[1:n2-2:2][1:n2-2:2] + tmp[3:n2-2:2][1:n2-2:2] + tmp[1:n2-2:2][3:n2-2:2] );
	temp2[1:n2-2][1:n2-2] /= 16.0;

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
	real **var = v;
	var[0:n][0:n] = val;
}

void MultiGrid::calc_res_to_temp()
{
	real **tmp = temp, **vt = v, **ft = f;

	tmp[1:n-2][1:n-2] = ft[1:n-2][1:n-2] - ( 4*vt[1:n-2][1:n-2] - ( vt[2:n-2][1:n-2] + vt[0:n-2][1:n-2] + vt[1:n-2][2:n-2] + vt[1:n-2][0:n-2] ) )/h2;

	// Set residual value at boundary to zero
	// tmp[0][0:n] = 0;
	// tmp[0:n][0] = 0;
	// tmp[n-1][0:n] = 0;
	// tmp[0:n][n-1] = 0;

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

void MultiGrid::copy_temp_to_v()
{
	copy(temp, v);
}

void MultiGrid::copy_v_to_temp()
{
	copy(v, temp);
}

void MultiGrid::copy_f_to_temp()
{
	copy(f, temp);
}

void MultiGrid::add_temp_to_v()
{
	real **var = v, **tmp = temp;
	var[1:n-2][1:n-2] = var[1:n-2][1:n-2] + tmp[1:n-2][1:n-2];
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