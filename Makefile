# Compile basicIBM with Intel c++ compiler

CC=icc
CFlags=-c -Wall
CLibs=-larmadillo

all: basicIBM

basicIBM: main.o cell.o
	$(CC) main.o cell.o -o basicIBM

main: main.cpp
	$(CC) main.cpp $(CFlags) $(Clibs)

cell: cell.cpp
	$(CC) cell.cpp $(CFlags) $(Clibs)

clean:
	rm -rf *.o
