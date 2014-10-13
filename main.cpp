#include "multigridsolver.h"

#define GRID_DENSITY_LEVEL 9
#define NUM_GRIDS 6
#define ACCURACY 1e-4

int main (int argc, char *argv[]) 
{
	MultiGridSolver solver(GRID_DENSITY_LEVEL, NUM_GRIDS);
	solver.solve(ACCURACY);
	return 0;
}