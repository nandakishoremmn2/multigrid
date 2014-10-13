CFLAGS = -Ofast -Wno-write-strings
CC = g++
SRC = main.cpp multigrid.cpp multigridsolver.cpp
OBJ = $(SRC:.cpp = .o)

test: $(OBJ)
	$(CC) $(CFLAGS) -o all $(OBJ)

clean:
	rm -f core *.o
