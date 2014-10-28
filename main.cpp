#include "multigridsolver.h"

#include <iostream>
#include <cstdlib>
#include <csignal>

MultiGridSolver *solver_ptr;

int *get_input(int argc, char *argv[]);
void signalHandler( int signum );

int main (int argc, char *argv[]) 
{
	int *input = get_input(argc, argv);

	MultiGridSolver solver(input);
	solver_ptr = &solver;

	signal(SIGINT, signalHandler);

	solver.solve(1e-9);
	return 0;
}

void signalHandler( int signum )
{
	solver_ptr->save();
	exit(0);  
}

int *get_input(int argc, char *argv[])
{
	int *data = new int[7];
	if( argc >= 6 && argc <= 7 ) // If cmd args are present
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
				std::cout<<"Enter v1, v2";
				std::cin>>data[3]>>data[4];
				break;
			case 1:
				std::cout<<"Enter v1, v2, mu";
				std::cin>>data[3]>>data[4]>>data[5];
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