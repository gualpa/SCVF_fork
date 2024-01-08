
#include "allvars.h"
#include <stdio.h>
//double RadIncrement;
//int    NumRanWalk;
//double BoxSize;          
//double MaxRadiusSearch;  
//double ProxyGridSize;  
//double FracRadius;       
//double DeltaThreshold;   
//double DeltaSeed;        
//double OverlapTol;      
//int    FormatTracers;
//int    NumFiles;
//char   FileTracers[MAXCHAR];      
//char   FileVoids[MAXCHAR];        
//int    OMPcores;         
//int    RSDist;           
//double Redshift;         
//double OmegaMatter;
//double OmegaLambda;
//double Hubble;           
//int    GDist;           
//double FidOmegaMatter;        
//double FidOmegaLambda;        
//double FidHubble;        
//int    WriteProfiles;    
//double MinProfileDist;   
//double MaxProfileDist;   
//int    NumProfileBins;   
//char   PathProfiles[MAXCHAR];     
//double InnerShellVel;  
//double OuterShellVel;  
//double ScalePos;
//double ScaleVel;
//int    RunFlag;

//vector <double> StepTime;
//vector <string> StepName;

//int    NumTrac;
//double MeanNumTrac;
//double MeanSeparation;
//int    NumVoid;
//double LBox[3];
//FILE   *logfile;


//scale_struct Scale;
//varConfiguration VarConfig;
//vector <tracers> Tracer;
vector <voids> Void;

void print_varConfigurtion(varConfiguration VarConfigAux, logs LogAux){
    fprintf(stdout,"\n Print varConfigurtion  \n" );
    fprintf(stdout,"RadIncrement= (%f) \n",VarConfigAux.RadIncrement);
    fflush(stdout);
    fprintf(stdout,"NumRanWalk= (%i) \n",VarConfigAux.NumRanWalk);
    fflush(stdout);
    fprintf(stdout,"BoxSize= (%f) \n",VarConfigAux.BoxSize);
    fflush(stdout);
    fprintf(stdout,"MaxRadiusSearch= (%f) \n",VarConfigAux.MaxRadiusSearch);
    fflush(stdout);
    fprintf(stdout,"ProxyGridSize= (%f) \n",VarConfigAux.ProxyGridSize);
    fflush(stdout);
    fprintf(stdout,"FracRadius= (%f) \n",VarConfigAux.FracRadius);
    fflush(stdout);
    fprintf(stdout,"DeltaThreshold = (%f) \n",VarConfigAux.DeltaThreshold);
    fflush(stdout);
    fprintf(stdout,"DeltaSeed = (%f) \n",VarConfigAux.DeltaSeed);
    fflush(stdout);
    fprintf(stdout,"OverlapTol = (%f) \n",VarConfigAux.OverlapTol);
    fflush(stdout);
    fprintf(stdout,"FormatTracers = (%i) \n",VarConfigAux.FormatTracers);
    fflush(stdout);
    fprintf(stdout,"NumFiles = (%i) \n",VarConfigAux.NumFiles);
    fflush(stdout);
    fprintf(stdout,"FileTracers = (%s) \n",VarConfigAux.FileTracers);
    fflush(stdout);
    fprintf(stdout,"FileVoids = (%s) \n",VarConfigAux.FileVoids);
    fflush(stdout);
    fprintf(stdout,"OMPcores = (%i) \n",VarConfigAux.OMPcores);
    fflush(stdout);
    fprintf(stdout,"RSDist = (%i) \n",VarConfigAux.RSDist);
    fflush(stdout);
    fprintf(stdout,"Redshift = (%f) \n",VarConfigAux.Redshift);
    fflush(stdout);
    fprintf(stdout,"OmegaMatter = (%f) \n",VarConfigAux.OmegaMatter);
    fflush(stdout);
    fprintf(stdout,"OmegaLambda = (%f) \n",VarConfigAux.OmegaLambda);
    fflush(stdout);
    fprintf(stdout,"Hubble = (%f) \n",VarConfigAux.Hubble);
    fflush(stdout);
    fprintf(stdout,"GDist = (%i) \n",VarConfigAux.GDist);
    fflush(stdout);
    fprintf(stdout,"FidOmegaMatter = (%f) \n",VarConfigAux.FidOmegaMatter);
    fflush(stdout);
    fprintf(stdout,"FidOmegaLambda = (%f) \n",VarConfigAux.FidOmegaLambda);
    fflush(stdout);
    fprintf(stdout,"FidHubble = (%f) \n",VarConfigAux.FidHubble);
    fflush(stdout);
    fprintf(stdout,"WriteProfiles = (%i) \n",VarConfigAux.WriteProfiles);
    fflush(stdout);
    fprintf(stdout,"MinProfileDist = (%f) \n",VarConfigAux.MinProfileDist);
    fflush(stdout);
    fprintf(stdout,"MaxProfileDist = (%f) \n",VarConfigAux.MaxProfileDist);
    fflush(stdout);
    fprintf(stdout,"NumProfileBins = (%i) \n",VarConfigAux.NumProfileBins);
    fflush(stdout);
    fprintf(stdout,"PathProfiles = (%s) \n",VarConfigAux.PathProfiles);
    fflush(stdout);
    fprintf(stdout,"InnerShellVel = (%f) \n",VarConfigAux.InnerShellVel);
    fflush(stdout);
    fprintf(stdout,"OuterShellVel = (%f) \n",VarConfigAux.OuterShellVel);
    fflush(stdout);
    fprintf(stdout,"RunFlag = (%i) \n",VarConfigAux.RunFlag);
    fflush(stdout);
    fprintf(stdout,"NumTrac = (%i) \n",VarConfigAux.NumTrac);
    fflush(stdout);
    fprintf(stdout,"MeanNumTrac = (%f) \n",VarConfigAux.MeanNumTrac);
    fflush(stdout);
    fprintf(stdout,"MeanSeparation = (%f) \n",VarConfigAux.MeanSeparation);
    fflush(stdout);
    fprintf(stdout,"NumVoid = (%f) \n",VarConfigAux.NumVoid);
    fflush(stdout);
    fprintf(stdout,"LBox = (%f) (%f) (%f)\n",VarConfigAux.LBox[0],VarConfigAux.LBox[1],VarConfigAux.LBox[2]);
    fflush(stdout);
    fprintf(stdout,"logfile = (%s) \n",LogAux.logfile);
    fflush(stdout);
    fprintf(stdout,"ScalePos Pos= (%f) \n",VarConfigAux.Scale.Pos);
    fflush(stdout);
    fprintf(stdout,"ScalePos Vel= (%f) \n",VarConfigAux.Scale.Vel);
    fflush(stdout);



} 
  