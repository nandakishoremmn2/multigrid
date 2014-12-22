#include <fstream>
#include "multigrid.h"
#include <cmath>
#include <cstdlib>

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h> 
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>

// GLOBALS
__global__ void d_set(real *v, real val, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x ;
    int j = blockIdx.y * blockDim.y + threadIdx.y ;

    if(i>=0 && i<n && j>=0 && j<n) {
    	v[i + n*j] = val;
    }
}

// assign a 2D distribution of CUDA "threads" within each CUDA "block"    


/*
	Constructors and destructors
*/
MultiGrid::MultiGrid(int grid_density, int no_of_child_grids, float wu, int lvl)
{
	n = pow(2, grid_density) + 1;

	WU = wu;
	level = lvl;

	h = pow(.5, grid_density);
	h2 = h*h;

	// Allocate memory for the grid on GPU
	v = initialise();
	f = initialise();
	temp = initialise();

	// Initialise the coarser child grids
	grid2 = ( no_of_child_grids > 0 ) ? new MultiGrid(grid_density-1, no_of_child_grids-1, WU/4.0, level+1) : NULL;

	ThreadsPerBlock=16;
	dimBlock = dim3( ThreadsPerBlock, ThreadsPerBlock );
	dimBlock1d = dim3( ThreadsPerBlock );

	// calculate number of blocks along X and Y in a 2D CUDA "grid"
	dimGrid = dim3( ceil(float(n)/float(dimBlock.x)), ceil(float(n)/float(dimBlock.y)), 1 );
	dimGrid1d = dim3( ceil(float(n)/float(dimBlock.x)), 1, 1 );


	apply_boundary_conditions();
	d_set<<<dimGrid, dimBlock>>>(v, 0., n);
	d_set<<<dimGrid, dimBlock>>>(f, 0., n);
	d_set<<<dimGrid, dimBlock>>>(temp, 0., n);

}

MultiGrid::~MultiGrid()
{
	// Free memory of variables on CPU
	deallocate(v);
	deallocate(f);
	deallocate(temp);

	// Destroy child grids
	delete grid2;
}

/*
	Private methods
*/


real *MultiGrid::initialise()
{
	real *d_var;
	cudaMalloc((void **)&d_var,n*n*sizeof(real));
	d_set<<<dimGrid, dimBlock>>>(d_var, 0., n);
	return d_var;
}

void MultiGrid::deallocate(real *d_var)
{
	cudaFree(d_var);
}

void MultiGrid::copy(real *src, real *dest)
{
	cudaMemcpy(dest, src, n*n*sizeof(real), cudaMemcpyDeviceToDevice);
}

real MultiGrid::norm2(real *data)
{
	// L-squared norm
	thrust::device_ptr<real> dptr(data);
	real norm2val = thrust::reduce(dptr, dptr + n*n, (real) 0, thrust::plus<real>());

	return sqrt(norm2val);
}

/*
	Public methods
*/

void MultiGrid::relax_once()
{
}
__global__ void d_relax_once_rb(real *v, real *f, int n, int nx, int ny, real h2)
{
    int i = 2 * ( blockIdx.x * blockDim.x + threadIdx.x ) + nx ;
    int j = 2 * ( blockIdx.y * blockDim.y + threadIdx.y ) + ny ;

    if(i>0 && i<n-1 && j>0 && j<n-1) {
		v[i + n*j] = ( v[i + n*(j+1)] + v[i + n*(j-1)] + v[i+1 + n*j] + v[i-1 + n*j] + h2*f[i + n*j])/4.0;
    }
}

void MultiGrid::relax_once_rb()
{
	d_relax_once_rb<<<dimGrid, dimBlock>>>(v, f, n, 0, 0, h2);
	d_relax_once_rb<<<dimGrid, dimBlock>>>(v, f, n, 0, 1, h2);
	d_relax_once_rb<<<dimGrid, dimBlock>>>(v, f, n, 1, 0, h2);
	d_relax_once_rb<<<dimGrid, dimBlock>>>(v, f, n, 1, 1, h2);
}

GridData MultiGrid::relax(int vn)
{
	for (int k = 0; k < vn; ++k)
	{
		relax_once_rb();
	}

	calc_res_to_temp();

	GridData data = {
		norm2(temp),
		vn*WU,
		vn,
		get_level()
		// 0,0,0,0
	};

	return data;
}

int MultiGrid::getSize()
{
	return n;
}

__global__ void d_interp(real *temp, real *temp2, int n, int n2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x ;
    int j = blockIdx.y * blockDim.y + threadIdx.y ;

    if(i>0 && i<n-1 && j>0 && j<n-1) {
		if (i%2==0)
		{
			if (j%2==0)
			{
				temp[i + n*j] = temp2[i/2 + n2*j/2];
			}
			else
			{
				temp[i + n*j] = ( temp2[i/2 + n2*(j-1)/2] + temp2[i/2 + n2*(j+1)/2] ) / 2.;
			}
		}
		else
		{
			if (j%2==0)
			{
				temp[i + n*j] = ( temp2[(i-1)/2 + n2*j/2] + temp2[(i+1)/2 + n2*j/2] ) / 2.;
			}
			else
			{
				temp[i + n*j] = ( temp2[(i-1)/2 + n2*(j-1)/2] + temp2[(i+1)/2 + n2*(j-1)/2] + temp2[(i-1)/2 + n2*(j+1)/2] + temp2[(i+1)/2 + n2*(j+1)/2] ) / 4.;
			}
		}
    }
}

void MultiGrid::interp()
{
	// Interpolate values from *grid2.temp to temp
	d_interp<<<dimGrid, dimBlock>>>(temp, grid2->temp, n, grid2->getSize());
}

__global__ void d_restrict(real *temp, real *temp2, int n, int n2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x ;
    int j = blockIdx.y * blockDim.y + threadIdx.y ;

    if(i>0 && i<n2-1 && j>0 && j<n2-1) {
		// Weighted average from all 9 points
		temp2[i + n2*j] = 4.0*temp[2*i + n*2*j];
		temp2[i + n2*j] += 2.0*( temp[2*i+1 + n*2*j] + temp[2*i-1 + n*2*j] + temp[2*i + n*(2*j-1)] + temp[2*i + n*(2*j+1)] );
		temp2[i + n2*j] += 1.0*( temp[2*i+1 + n*(2*j+1)] + temp[2*i-1 + n*(2*j-1)] + temp[2*i+1 + n*(2*j-1)] + temp[2*i-1 + n*(2*j+1)] );
		temp2[i + n2*j] /= 16.0;
    }
}

void MultiGrid::restrict()
{
	// Restrict values from temp to *grid2.temp
	d_restrict<<<grid2->dimGrid, grid2->dimBlock>>>(temp, grid2->temp, n, grid2->getSize());

}

__global__ void d_apply_boundary_conditions(real *v, real h2, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x ;

    if(i>=0 && i<n) {
    	v[0 + n*i] = v[i + 0] = (real)i*i*h2;
		v[n-1 + n*i] = v[i + n*(n-1)] = 1.0 - (real)i*i*h2;
    }
}

void MultiGrid::apply_boundary_conditions()
{
	d_apply_boundary_conditions<<<dimGrid1d, dimBlock1d>>>(v, h2, n);
}

void MultiGrid::set_v(real val)
{
	d_set<<<dimGrid, dimBlock>>>(v, val, n);
}

__global__ void d_calc_res_to_temp(real *temp, real *v, real *f, int n, real h2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x ;
    int j = blockIdx.y * blockDim.y + threadIdx.y ;

    if(i>0 && i<n-1 && j>0 && j<n-1) {
		temp[i + n*j] = f[i + n*j] - ( 4*v[i + n*j] - ( v[i+1 + n*j] + v[i-1 + n*j] + v[i + n*(j+1)] + v[i + n*(j-1)] ) )/h2;
    }
    else if (i==0 || j==0 || i==n-1 || j==n-1)
    {
    	temp[i + n*j] = 0;
    }
}

void MultiGrid::calc_res_to_temp()
{
	d_calc_res_to_temp<<<dimGrid, dimBlock>>>(temp, v, f, n, h2);
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

__global__ void d_add_temp_to_v(real *temp, real *v, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x ;
    int j = blockIdx.y * blockDim.y + threadIdx.y ;

    if(i>0 && i<n-1 && j>0 && j<n-1) {
		v[i + n*j] = v[i + n*j] + temp[i + n*j];
    }
}

void MultiGrid::add_temp_to_v()
{
	d_add_temp_to_v<<<dimGrid, dimBlock>>>(temp, v, n);
}

void MultiGrid::save_grid(char *filename)
{
	std::ofstream outfile(filename);
	real *temp_v = new real[n*n];
	cudaMemcpy(temp_v, v, n*n*sizeof(real), cudaMemcpyDeviceToHost);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			outfile<<temp_v[i + n*j]<<" ";
		}
		outfile<<"\n";
	}
	delete [] temp_v;
}

real MultiGrid::get_L2norm()
{
	calc_res_to_temp();
	return norm2(temp);
}

real MultiGrid::get_wu()
{
	return WU;
}

int MultiGrid::get_level()
{
	return level;
}