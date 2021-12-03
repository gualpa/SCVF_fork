#__________________________________________________________
#                                                  Compiler
                                                  
CC = g++

#__________________________________________________________
#                                                    voro++

VORO_LIB = -L/usr/local/lib -lvoro++
VORO_INC = -I/usr/local/include/voro++

#__________________________________________________________
#                                                    OpenMP

OMP = -fopenmp

#__________________________________________________________
#

EXEC = main.x

OBJS = allvars.o tools.o qromb.o io.o voronoi.o grid.o finder.o \
       qsort.o velocity.o profiles.o main.o

INCL = allvars.h proto.h Makefile

CFLAGS += $(OMP) $(VORO_INC) 
   	
LIBS = $(VORO_LIB) -lm

$(EXEC): $(OBJS) 
	$(CC) $(OMP) $(OBJS) $(LIBS) -o $(EXEC)  

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS) $(EXEC) *~ 
