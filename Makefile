#__________________________________________________________
#                                                  Compiler
                                                  
CPP = g++

#__________________________________________________________
#                                                    voro++

VORO_LIB = -L/home/srgualpa/projects/SCVF_Polaco_finder_void/voro++-0.4.6/src -lvoro++
VORO_INC = -I/home/srgualpa/projects/SCVF_Polaco_finder_void/voro++-0.4.6/src

#__________________________________________________________
#                                                       GSL

GSL_LIB = -L/opt/spack/dev/opt/spack/linux-centos7-x86_64_v3/gcc-12.1.0/gsl-2.7.1-ttp5l7ez636aw4kkpbfzfufdu2546jgf/lib -lgsl -lgslcblas 
GSL_INC = -I/opt/spack/dev/opt/spack/linux-centos7-x86_64_v3/gcc-12.1.0/gsl-2.7.1-ttp5l7ez636aw4kkpbfzfufdu2546jgf/include/gsl/


#__________________________________________________________
#                                                    OpenMP

OMP = -fopenmp

#__________________________________________________________
#

EXEC = main.x

OBJS = allvars.o tools.o cosmology.o io.o voronoi.o grid.o finder.o \
       qsort.o velocity.o profiles.o main.o

INCL = allvars.h tools.h cosmology.h io.h voronoi.h grid.h finder.h \
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
