# Compile basicIBM with Intel c++ compiler

CC=icpc
CFlags=-c -Wall -std=c++11 -ggdb -O0
CLibs=-larmadillo
LDFLAGS=-ggdb -std=c++11 -O0 -Wall

all: basicIBM clean

basicIBM: main.o element.o cell.o tissue.o BoundToGrid1.o BoundToGrid2.o GridToBound.o NavierStokes.o GlobalToLocal.o smallfunctions.o ReadParams.o
	$(CC) main.o element.o cell.o tissue.o BoundToGrid1.o BoundToGrid2.o GridToBound.o NavierStokes.o GlobalToLocal.o smallfunctions.o ReadParams.o -o basicIBM $(LDFLAGS)

main: main.cpp
	$(CC) $(CFlags) main.cpp $(CLibs)

element: element.cpp
	$(CC) $(CFlags) element.cpp $(CLibs)

cell: cell.cpp
	$(CC) $(CFlags) cell.cpp $(CLibs)

tissue: tissue.cpp
	$(CC) $(CFlags) tissue.cpp $(CLibs)

BoundToGrid1: BoundToGrid1.cpp
	$(CC) $(CFlags) BoundToGrid1.cpp $(CLibs)

BoundToGrid2: BoundToGrid2.cpp
	$(CC) $(CFlags) BoundToGrid2.cpp $(CLibs)

GridToBound: GridToBound.cpp
	$(CC) $(CFlags) GridToBound.cpp $(CLibs)

NavierStokes: NavierStokes.cpp
	$(CC) $(CFlags) NavierStokes.cpp $(CLibs)

GlobalToLocal: GlobalToLocal.cpp
	$(CC) $(CFlags) GlobalToLocal.cpp $(CLibs)

ReadParams: ReadParams.cpp
	$(CC) $(CFlags) ReadParams.cpp $(CLibs)

smallfunctions: smallfunctions.cpp
	$(CC) $(CFlags) smallfunctions.cpp $(CLibs)

clean:
	rm -rf *.o
