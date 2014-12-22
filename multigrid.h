#ifndef MULTIGRID_H
#define MULTIGRID_H

#include <vector_types.h>

typedef double real;

/**
	This structure is to be returned by the relax()
*/
struct GridData
{
	real residue;
	real wu;
	int iterations;
	int level;
};

class MultiGrid
{
public:

	/**
		Initialises the multigrid
		@param grid_density The no. of grid points = ( 2^grid_density + 1 )^2
		@param no_of_child_grids The no. of coarser grids it has ( < grid_density )
		@param wu Work unit cost for one iteration
		@param lvl The grid level ( 0 for coarsest )
	*/
	MultiGrid(int grid_density, int no_of_child_grids, float wu, int lvl); 

	~MultiGrid();

	/**
		Iterates the grid a given no. of times
		@param vn the no. of iterations to perform
		@return returns the norm-2 of residue
	*/
	GridData relax(int vn);

	/**
		Relaes the grid once
	*/
	void relax_once();

	/**
		Relaes the grid once (red-black)
	*/
	void relax_once_rb();

	/**
		Returns the size of the grid
		@return The size of the grid
	*/
	int getSize();

	/**
		Contains the pointer to the next coarser grid
		The value is NULL if this grid is the coarsest grid
	*/
	MultiGrid *grid2;

	/**
		Applies interpolation from the coarser grid (i.e. grid2) to this grid
		More specifically it transfers values between the temporary variables only
	*/
	void interp();

	/**
		Applies restriction from this grid to the coarser grid (i.e. grid2)
		More specifically it transfers values between the temporary variables only
	*/
	void restrict();

	/**
		Assigns the boundary values
	*/
	void apply_boundary_conditions();

	/**
		Sets the v to the given value
		@param val The value to set
	*/
	void set_v(real);

	/**
		Calculates the residue and copies it to temp variable
	*/
	void calc_res_to_temp();

	/**
		Copies temp values to f
	*/
	void copy_temp_to_f();

	/**
		Copies temp values to v
	*/
	void copy_temp_to_v();

	/**
		Copies v values to temp
	*/
	void copy_v_to_temp();

	/**
		Copies f values to temp
	*/
	void copy_f_to_temp();

	/**
		Adds the temp values to existing v values ( v = v + temp )
	*/
	void add_temp_to_v();

	/**
		Saves the grid data to a file
	*/
	void save_grid(char *filename);

	/**
		Temporary variable
	*/
	real *temp;
	real *d_temp;

	/**
		Returns the L squared norm of residue
			Calculate redidue and stores to temp and
			then finds the norm2 of temp
	*/
	real get_L2norm();

	/**
		Gets work units per iteration of the grid
	*/
	real get_wu();

	/**
		Gets the grid level
	*/
	int get_level();

	int ThreadsPerBlock;
	dim3 dimBlock;
	dim3 dimBlock1d;

	// calculate number of blocks along X and Y in a 2D CUDA "grid"
	dim3 dimGrid;
	dim3 dimGrid1d;

// private:
	/**
		The size of the grid
	*/
	int n;

	/** 
		Grid spacing 
	*/
	real h;
	/** 
		Square of h 
	*/
	real h2;

	/**
		Stores the solution at all levels
	*/
	real *v;

	/**
		RHS of the equation Au=f
	*/
	real *f;

	/**
		Work unit cost for each itertion
	*/
	float WU;

	/**
		Grid level
	*/
	float level;

	/**
		Initialises and allocates memory for the pointer variable passed
		@param var pointer to the variable
	*/
	real *initialise();

	/**
		Does the opposite of initialise
		@param var pointer to freed
		@see initialise
	*/
	void deallocate(real *var);

	/**
		Copies values from one grid to another
		@param src the source grid
		@param dest the destination grid
	*/
	void copy(real *src, real *dest);

	/**
		Calculates norm-2 of the data provided
		@param data the source data
		@return norm2 the norm-2 of the data
	*/
	real norm2(real *data);
};

#endif // MULTIGRID_H