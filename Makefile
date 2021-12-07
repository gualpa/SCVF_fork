#__________________________________________________________
#                                                  Compiler
                                                  
CPP = g++

#__________________________________________________________
#                                                    voro++

VORO_LIB = -L/usr/local/lib -lvoro++
VORO_INC = -I/usr/local/include/voro++

#__________________________________________________________
#                                                       GSL

GSL_LIB = -L/usr/lib -lgsl -lgslcblas 
GSL_INC = -I/usr/include/gsl

#__________________________________________________________
#                                                    OpenMP

OMP = -fopenmp

#__________________________________________________________
#

EXEC = main.x

OBJS = allvars.o tools.o qromb.o io.o voronoi.o grid.o finder.o \
       qsort.o velocity.o profiles.o main.o

INCL = allvars.h tools.h qromb.h io.h voronoi.h grid.h finder.h \
       qsort.h velocity.h profiles.h Makefile

CFLAGS = $(OMP) $(VORO_INC) $(GSL_INC) 
   	
LIBS = $(VORO_LIB) $(GSL_LIB) -lm

build: $(EXEC)

$(EXEC): $(OBJS) 
	$(CPP) $(OMP) $(OBJS) $(LIBS) -o $(EXEC)  

%.o: %.cc $(INCL)
	$(CPP) $(CFLAGS) -c $< -o $@

#$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS) $(EXEC) *~ 
