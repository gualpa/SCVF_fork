/**
 * @file allvars.h
 * @brief Archivo con definición tipos de datos especiales para el programa
 *
 * Copyright 2023 Andrés Nicolás Ruiz, Sebastián Rogelio Gualpa
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may
 * be used to endorse or promote products derived from this software without specific
 * prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
// allvars.h
#ifndef ALLVARS_H
#define ALLVARS_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <vector>
#include "voro++.hh"

using namespace std;
using namespace voro;

#define MAXCHAR 500

// Some constants

#define PI   3.141592653589793238
#define CVEL  299792.469 

struct scale_struct {
  double Pos;
  double Vel;
};

struct varConfiguration {
  double RadIncrement;
  int    NumRanWalk;
  double BoxSize;          
  double MaxRadiusSearch;  
  double ProxyGridSize;  
  double FracRadius;       
  double DeltaThreshold;   
  double DeltaSeed;        
  double OverlapTol;   
  int    FormatTracers;
  int    NumFiles;
  char   FileTracers[MAXCHAR];      
  char   FileVoids[MAXCHAR];        
  int    OMPcores;         
  int    RSDist;           
  double Redshift;         
  double OmegaMatter;
  double OmegaLambda;
  double Hubble;           
  int    GDist;           
  double FidOmegaMatter;        
  double FidOmegaLambda;        
  double FidHubble;        
  int    WriteProfiles;    
  double MinProfileDist;   
  double MaxProfileDist;   
  int    NumProfileBins;   
  char   PathProfiles[MAXCHAR];     
  double InnerShellVel;  
  double OuterShellVel;  
  double ScalePos;
  double ScaleVel;
  int    RunFlag;
  int    NumTrac;
  double MeanNumTrac;
  double MeanSeparation;
  int    NumVoid;
  double LBox[3];
  scale_struct Scale;
}; 

// Tracers

struct tracers {
  float Pos[3];
  float Vel[3];
  float Cen[3];
  float Delta;
  float Volume;
}; 

// Voids

struct voids {
  float Rad;
  float Rini;
  float Ini[3];
  float Pos[3];
  float Vel[3];
  float Dtype;
  float Delta;
  float Poisson;
  bool  ToF;  
  int   Nran;
};

struct  logs{
  FILE   *logfile;
  vector <double> StepTime;
  vector <string> StepName;
};

void print_varConfigurtion(varConfiguration VarConfigAux, logs LogAux);


#endif // ALLVARS_H