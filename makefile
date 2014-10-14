CFLAGS = -Ofast -Wno-write-strings
CILKFLAGS = -fcilkplus
OPENMPFLAGS = -fopenmp
CC = g++-4.9
SRC = main.cpp multigrid.cpp multigridsolver.cpp
OBJ = $(SRC:.cpp = .o)

all: $(OBJ)
	$(CC) $(CFLAGS) -o solver_serial $(OBJ)

cilk: $(OBJ)
	$(CC) $(CFLAGS) $(CILKFLAGS) -o solver_cilk $(OBJ)

openmp: $(OBJ)
	$(CC) $(CFLAGS) $(OPENMPFLAGS) -o solver_openmp $(OBJ)

clean:
	rm -f core *.o