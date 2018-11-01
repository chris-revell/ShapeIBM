# Compile basicIBM with Intel c++ compiler

CC=icc
CFlags=-c -Wall -std=c++11 -g -O0
CLibs=-larmadillo

all: basicIBM clean

basicIBM: main.o cell.o tissue.o BoundToGrid1.o BoundToGrid2.o GridToBound.o NavierStokes.o GlobalToLocal.o smallfunctions.o ReadParams.o
	$(CC) main.o cell.o tissue.o BoundToGrid1.o BoundToGrid2.o GridToBound.o NavierStokes.o GlobalToLocal.o smallfunctions.o ReadParams.o -o basicIBM

main: main.cpp
	$(CC) main.cpp $(CFlags) $(CLibs)

cell: cell.cpp
	$(CC) cell.cpp $(CFlags) $(CLibs)

tissue: tissue.cpp
	$(CC) tissue.cpp $(CFlags) $(CLibs)

BoundToGrid1: BoundToGrid1.cpp
	$(CC) BoundToGrid1.cpp $(CFlags) $(CLibs)

BoundToGrid2: BoundToGrid2.cpp
	$(CC) BoundToGrid2.cpp $(CFlags) $(CLibs)

GridToBound: GridToBound.cpp
	$(CC) GridToBound.cpp $(CFlags) $(CLibs)

NavierStokes: NavierStokes.cpp
	$(CC) NavierStokes.cpp $(CFlags) $(CLibs)

GlobalToLocal: GlobalToLocal.cpp
	$(CC) GlobalToLocal.cpp $(CFlags) $(CLibs)

smallfunctions: smallfunctions.cpp
	$(CC) smallfunctions.cpp $(CFlags) $(CLibs)

	ReadParams: ReadParams.cpp
		$(CC) ReadParams.cpp $(CFlags) $(CLibs)

clean:
	rm -rf *.o
