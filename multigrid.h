#ifndef MULTIGRID_H
#define MULTIGRID_H

typedef double real;

class MultiGrid
{
public:

	/**
		Initialises the multigrid
		@param grid_density The no. of grid points = ( 2^grid_density + 1 )^2
		@param no_of_child_grids The no. of coarser grids it has ( < grid_density )
	*/
	MultiGrid(int grid_density, int no_of_child_grids); 

	~MultiGrid();

	/**
		Iterates the grid a given no. of times
		@param v the no. of iterations to perform
	*/
	void relax(int v);

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
		Calculates the residue and copies it to temp variable
	*/
	void calc_res_to_temp();

	/**
		Copies temp values to f
	*/
	void copy_temp_to_f();

	/**
		Copies v values to temp
	*/
	void copy_v_to_temp();

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
	real **temp;

private:
	/**
		The size of the grid
	*/
	int n;

	/**
		Stores the solution at all levels
	*/
	real **v;

	/**
		RHS of the equation Au=f
	*/
	real **f;

	/**
		Copies values from one grid to another
		@param src the source grid
		@param dest the destination grid
	*/
	void copy(real **src, real **dest);

};

#endif // MULTIGRID_H