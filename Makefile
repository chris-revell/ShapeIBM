# Compile basicIBM with Intel c++ compiler

#CC=icpc
#CFlags=-c -Wall -std=c++11 -ggdb -O0
#CLibs=-larmadillo
#LDFLAGS=-ggdb -std=c++11 -O0 -Wall

CC := icpc
SRCDIR := src
BUILDDIR := build
TARGET := basicIBM

#MKLROOT := /opt/intel/compilers_and_libraries_2018.3.185/mac/compiler
MKLROOT := /opt/intel/mkl
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -c -Wall -std=c++11 -ggdb -O0 -fopenmp -mkl=parallel  # -Wall
LIB := -larmadillo -L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/../compiler/lib -mkl#-liomp5 -lpthread -lm -ldl -L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/../compiler/lib -mkl
INC := -I include

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

#clean:
#  @echo " Cleaning...";
#  @echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)
#
#
#all: basicIBM clean
#
#basicIBM: main.o element.o cell.o tissue.o BoundToGrid1.o BoundToGrid2.o GridToBound.o NavierStokes.o GlobalToLocal.o smallfunctions.o ReadParams.o
#	$(CC) main.o element.o cell.o tissue.o BoundToGrid1.o BoundToGrid2.o GridToBound.o NavierStokes.o GlobalToLocal.o smallfunctions.o ReadParams.o -o basicIBM $(LDFLAGS)
#
#main: main.cpp
#	$(CC) $(CFlags) main.cpp $(CLibs)
#
#element: element.cpp
#	$(CC) $(CFlags) element.cpp $(CLibs)
#
#cell: cell.cpp
#	$(CC) $(CFlags) cell.cpp $(CLibs)
#
#tissue: tissue.cpp
#	$(CC) $(CFlags) tissue.cpp $(CLibs)
#
#BoundToGrid1: BoundToGrid1.cpp
#	$(CC) $(CFlags) BoundToGrid1.cpp $(CLibs)
#
#BoundToGrid2: BoundToGrid2.cpp
#	$(CC) $(CFlags) BoundToGrid2.cpp $(CLibs)
#
#GridToBound: GridToBound.cpp
#	$(CC) $(CFlags) GridToBound.cpp $(CLibs)
#
#NavierStokes: NavierStokes.cpp
#	$(CC) $(CFlags) NavierStokes.cpp $(CLibs)
#
#GlobalToLocal: GlobalToLocal.cpp
#	$(CC) $(CFlags) GlobalToLocal.cpp $(CLibs)
#
#ReadParams: ReadParams.cpp
#	$(CC) $(CFlags) ReadParams.cpp $(CLibs)
#
#smallfunctions: smallfunctions.cpp
#	$(CC) $(CFlags) smallfunctions.cpp $(CLibs)
#
#clean:
#	rm -rf *.o
#
