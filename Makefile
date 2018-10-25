# Compile basicIBM with Intel c++ compiler

CC=icc
CFlags=-c -Wall -std=c++11
CLibs=-larmadillo

all: basicIBM clean

basicIBM: main.o cell.o
	$(CC) main.o cell.o -o basicIBM

main: main.cpp
	$(CC) main.cpp $(CFlags) $(CLibs)

cell: cell.cpp
	$(CC) cell.cpp $(CFlags) $(CLibs)

clean:
	rm -rf *.o
