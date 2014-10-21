#include "multigridsolver.h"

#include <iostream>
#include <cstdlib>

#define GRID_DENSITY_LEVEL 9
#define NUM_GRIDS 6
#define ACCURACY 1e-4

int *get_input(int argc, char *argv[]);

int main (int argc, char *argv[]) 
{
	int *input = get_input(argc, argv);

	MultiGridSolver solver(input);
	solver.solve(ACCURACY);
	return 0;
}

int *get_input(int argc, char *argv[])
{
	int *data = new int[7];
	if( argc >= 7 && argc <= 8 ) // If cmd args are present
	{
		for (int i = 0; i < argc-1; ++i)
		{
			data[i] = atoi(argv[i+1]);
		}
	}
	else 
	{
		std::cout<<"Enter Grid density (no. of points = (2^(grid density)+1)^2 :";
		std::cin>>data[0];
		std::cout<<"Enter no. of grids to use ( < grid density) :";
		std::cin>>data[1];	
		std::cout<<"Enter scheme (0, V cycle) (1, Mu cycle), (2, FMG) : ";
		std::cin>>data[2];
		switch(data[2])
		{
			case 0:
				std::cout<<"Enter v1, v2, cycles : ";
				std::cin>>data[3]>>data[4]>>data[5];
				break;
			case 1:
				std::cout<<"Enter v1, v2, mu, cycles : ";
				std::cin>>data[3]>>data[4]>>data[5]>>data[6];
				break;
			case 2:
				std::cout<<"Enter v0, v1, v2 : ";
				std::cin>>data[3]>>data[4]>>data[5];
				break;
			// default:
				// show_help_and_exit();
		}
	}

	return data;
}