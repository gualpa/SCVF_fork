#
# @file allvars.h
# @brief Archivo Makefile para compilar y ejecutar el programa que identifica voids astronómicos.
#
# Copyright 2023 Andrés Nicolás Ruiz, Sebastián Rogelio Gualpa
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may
# be used to endorse or promote products derived from this software without specific
# prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#

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
