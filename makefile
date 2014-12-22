CFLAGS = -O3 -Wno-write-strings
CILKFLAGS = -fcilkplus
OPENMPFLAGS = -fopenmp
CC = g++
SRC = main.cpp multigrid.cpp multigridsolver.cpp
OBJ = $(SRC:.cpp = .o)

all: $(OBJ)
	$(CC) $(CFLAGS) -o solver_serial $(OBJ)

cilk: $(OBJ)
	$(CC) $(CFLAGS) $(CILKFLAGS) -o solver_cilk $(OBJ)

openmp: $(OBJ)
	$(CC) $(CFLAGS) $(OPENMPFLAGS) -o solver_openmp $(OBJ)

cuda:
	nvcc -O3 main.cpp multigridsolver.cpp multigrid.cu -o solver_cuda

clean:
	rm -f core *.o a.out out.dat data.dat nohup.out