
#include "allvars.h"
#include "io.h"
#include "tools.h"
#include "qromb.h"

void ReadInputFile(char *filename)
{
#define DOUBLE  1
#define STRING  2
#define INT     3
#define MAXTAGS 300

  int  i,j,nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  FILE *fd;
  char buf[MAXCHAR],buf1[MAXCHAR];
  char buf2[MAXCHAR],buf3[MAXCHAR];
  char fname[MAXCHAR];

  nt = 0;

  strcpy(tag[nt],"BoxSize");
  addr[nt] = &BoxSize;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"MaxRadiusSearch");
  addr[nt] = &MaxRadiusSearch;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"DeltaThreshold");
  addr[nt] = &DeltaThreshold;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"DeltaSeed");
  addr[nt] = &DeltaSeed;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"OverlapTol");
  addr[nt] = &OverlapTol;
  id[nt++] = DOUBLE;  

  strcpy(tag[nt],"FormatTracers");
  addr[nt] = &FormatTracers;
  id[nt++] = INT;

  strcpy(tag[nt],"NumFiles");
  addr[nt] = &NumFiles;
  id[nt++] = INT;

  strcpy(tag[nt],"FileTracers");
  addr[nt] = FileTracers;
  id[nt++] = STRING;

  strcpy(tag[nt],"FileVoids");
  addr[nt] = FileVoids;
  id[nt++] = STRING;

  strcpy(tag[nt],"ScalePos");
  addr[nt] = &ScalePos;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt],"ScaleVel");
  addr[nt] = &ScaleVel;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"OMPcores");
  addr[nt] = &OMPcores;
  id[nt++] = INT;

  strcpy(tag[nt],"RadIncrement");
  addr[nt] = &RadIncrement;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt],"NumRanWalk");
  addr[nt] = &NumRanWalk;
  id[nt++] = INT;

  strcpy(tag[nt],"FracRadius");
  addr[nt] = &FracRadius;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"RSDist");
  addr[nt] = &RSDist;
  id[nt++] = INT;

  strcpy(tag[nt],"Redshift");
  addr[nt] = &Redshift;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"OmegaMatter");
  addr[nt] = &OmegaMatter;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"OmegaLambda");
  addr[nt] = &OmegaLambda;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"Hubble");
  addr[nt] = &Hubble;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"GDist");
  addr[nt] = &GDist;
  id[nt++] = INT;

  strcpy(tag[nt],"FidOmegaMatter");
  addr[nt] = &FidOmegaMatter;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"FidOmegaLambda");
  addr[nt] = &FidOmegaLambda;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"FidHubble");
  addr[nt] = &FidHubble;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"WriteProfiles");
  addr[nt] = &WriteProfiles;
  id[nt++] = INT;

  strcpy(tag[nt],"MinProfileDist");
  addr[nt] = &MinProfileDist;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"MaxProfileDist");
  addr[nt] = &MaxProfileDist;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"NumProfileBins");
  addr[nt] = &NumProfileBins;
  id[nt++] = INT;

  strcpy(tag[nt],"PathProfiles");
  addr[nt] = PathProfiles;
  id[nt++] = STRING;

  strcpy(tag[nt],"InnerShellVel");
  addr[nt] = &InnerShellVel;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"OuterShellVel");
  addr[nt] = &OuterShellVel;
  id[nt++] = DOUBLE; 

  fd = SafeOpen(filename,"r");

  sprintf(fname,"%s.log",filename);	  
  logfile = SafeOpen(fname,"w");

  fprintf(logfile,"\n CONFIGURATION PARAMETERS USED \n\n");

  while (!feof(fd)) {

      buf[0] = 0;
      fgets(buf, 200, fd);

      if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2) continue;
      if (buf1[0] == '%') continue;

      for (i = 0, j = -1; i < nt; i++) {
          if (strcmp(buf1, tag[i]) == 0) {
    	     j = i;
    	     tag[i][0] = 0;
    	     break;
          }
      }

      if (j >= 0) {

         switch (id[j]) {

    	    case DOUBLE:
    	    *((double *) addr[j]) = atof(buf2);
    	    fprintf(logfile, " %-35s%g\n", buf1, *((double *) addr[j]));
    	    break;

    	    case STRING:
    	    strcpy((char *)addr[j], buf2);
    	    fprintf(logfile, " %-35s%s\n", buf1, buf2);
    	    break;

    	    case INT:
    	    *((int *) addr[j]) = atoi(buf2);
    	    fprintf(logfile, " %-35s%d\n", buf1, *((int *) addr[j]));
    	    break;

    	 }

      } else {

    	fprintf(stdout, "Error in file %s:  Tag '%s' not allowed or multiple defined. \n", filename, buf1);
        fflush(stdout);
	exit(EXIT_FAILURE);

      }
  }
  
  fclose(fd);

  for (i=0; i<nt; i++) {
      if (*tag[i]) {
	 fprintf(stdout, "Error. Missing value for tag '%s' in parameter file '%s'.\n", tag[i], filename);
         fflush(stdout);
	 exit(EXIT_FAILURE);
      }
  }

  fflush(logfile);

#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
}

void ReadTracers()
{

   int     i;
   FILE    *fd;
   double   Volume;
   clock_t t;

   fprintf(logfile,"\n READING TRACERS \n");
   t = clock();

   switch (FormatTracers) {

      case 0: 
      fprintf(logfile," | Reading ASCII format \n");	   
      ReadTracers_ASCII();
      break;

      case 1:
      fprintf(logfile," | Reading GADGET-1 format \n");	   
      ReadTracers_GADGET1();
      break;

      case 2:
      fprintf(logfile," | Reading GADGET-2 format \n");	   
      ReadTracers_GADGET2();
      break;

      case 3:
      fprintf(logfile," | Reading MXXL format \n");	   
      ReadTracers_MXXL();
      break;

      case 4:
      fprintf(logfile," | Reading BINARY format \n");	   
      ReadTracers_BINARY();

   }

   LBox[0] = BoxSize;
   LBox[1] = BoxSize;
   LBox[2] = BoxSize;

   if (RSDist == 1) RedshiftSpaceDistortions();
   if (GDist == 1) GeometricalDistortions();

   Volume = LBox[0]*LBox[1]*LBox[2];
   MeanNumTrac = (double)NumTrac/Volume;
   MeanSeparation = cbrt(Volume/(double)NumTrac);
 
   fprintf(logfile," | Number of tracers = %d \n",NumTrac);
   fprintf(logfile," | Size of the box: x-axis = %f \n",LBox[0]);
   fprintf(logfile," |                  y-axis = %f \n",LBox[1]);
   fprintf(logfile," |                  z-axis = %f \n",LBox[2]);
   fprintf(logfile," | Mean number density [h³/Mpc³] = %e \n",MeanNumTrac);
   fprintf(logfile," | Mean separation [Mpc/h] = %e \n",MeanSeparation);

   StepName.push_back("Reading tracers");
   StepTime.push_back(Time(t,1));

}

void ReadTracers_ASCII()
{
   int  i;
   FILE *fd;

   NumTrac = CountLines(FileTracers);
   Tracer = (struct tracers *) malloc(NumTrac*sizeof(struct tracers));

   fd = SafeOpen(FileTracers,"r");
   for (i=0; i<NumTrac; i++) {
       fscanf(fd,"%f %f %f %f %f %f \n",&Tracer[i].Pos[0],&Tracer[i].Pos[1],&Tracer[i].Pos[2],	   
                                        &Tracer[i].Vel[0],&Tracer[i].Vel[1],&Tracer[i].Vel[2]);
       Tracer[i].Pos[0] *= ScalePos;
       Tracer[i].Pos[1] *= ScalePos;
       Tracer[i].Pos[2] *= ScalePos;
       Tracer[i].Vel[0] *= ScaleVel;
       Tracer[i].Vel[1] *= ScaleVel;
       Tracer[i].Vel[2] *= ScaleVel;
       //fscanf(fd,"%f %f %f \n",&Tracer[i].Pos[0],&Tracer[i].Pos[1],&Tracer[i].Pos[2]);	   
       //Tracer[i].Pos[0] *= ScalePos;
       //Tracer[i].Pos[1] *= ScalePos;
       //Tracer[i].Pos[2] *= ScalePos;
       //Tracer[i].Vel[0] = 0.0;
       //Tracer[i].Vel[1] = 0.0;
       //Tracer[i].Vel[2] = 0.0;
   }
   fclose(fd);

}

void ReadTracers_BINARY()
{
   int   i;
   FILE  *fd;

   fd = SafeOpen(FileTracers,"r");
   
   fread(&NumTrac,sizeof(int),1,fd);
   Tracer = (struct tracers *) malloc(NumTrac*sizeof(struct tracers));

   for (i=0; i<NumTrac; i++) {
       fread(&Tracer[i].Pos[0],sizeof(float),3,fd);
       fread(&Tracer[i].Vel[0],sizeof(float),3,fd);
       Tracer[i].Pos[0] *= ScalePos;
       Tracer[i].Pos[1] *= ScalePos;
       Tracer[i].Pos[2] *= ScalePos;
       Tracer[i].Vel[0] *= ScaleVel;
       Tracer[i].Vel[1] *= ScaleVel;
       Tracer[i].Vel[2] *= ScaleVel;
   }
   fclose(fd);

}

void ReadTracers_GADGET1()
{

  struct GadgetHeader {
    int      Npart[6];
    double   Mass[6];
    double   Time;
    double   Redshift;
    int      Flag_1;
    int      Flag_2;
    int      NpartTotal[6];
    int      Flag_3;
    int      NumFiles;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam; 
    char     fill[96];  
  } Header;

  int   i,j,k,SkipBlock,dummy,Np,id,NC;
  float pos[3],vel[3];
  FILE  *f1,*f2,*f3;
  char  snapshot[MAXCHAR];

  if (NumFiles == 1) 
     sprintf(snapshot,"%s",FileTracers);	  
  else 
     sprintf(snapshot,"%s.0",FileTracers);	  

  f1 = SafeOpen(snapshot,"r");

  fread(&dummy,sizeof(int),1,f1);
  fread(&Header,sizeof(struct GadgetHeader),1,f1); 
  fclose(f1);

  NumTrac = Header.NpartTotal[1];
  
  if (Header.BoxSize*ScalePos != BoxSize || Header.NumFiles != NumFiles) {
     fprintf(stdout,"\nError. Missmatch with Gadget header.\n");
     fprintf(stdout,"BoxSize = %f (%f in inputfile)\n",Header.BoxSize*ScalePos,BoxSize);
     fprintf(stdout,"NumFiles = %d (%d in inputfile)\n",Header.NumFiles,NumFiles);
     fflush(stdout);
     exit(EXIT_FAILURE);
  }

  if ((RSDist == 1 || GDist == 1) && (Header.Omega0 != OmegaMatter || 
       Header.OmegaLambda != OmegaLambda || Header.Redshift != Redshift || 
       Header.HubbleParam*100.0 != Hubble)) { 
     fprintf(stdout,"\nError. Missmatch with Gadget header.\n");
     fprintf(stdout,"OmegaMatter = %f (%f in inputfile)\n",Header.Omega0,OmegaMatter);
     fprintf(stdout,"OmegaLambda = %f (%f in inputfile)\n",Header.OmegaLambda,OmegaLambda);
     fprintf(stdout,"Hubble = %f (%f in inputfile)\n",Header.HubbleParam*100.0,Hubble);
     fprintf(stdout,"Redshift = %f (%f in inputfile)\n",Header.Redshift,Redshift);
     fflush(stdout);
     exit(EXIT_FAILURE);
  }
  
  Tracer = (struct tracers *) malloc(NumTrac*sizeof(struct tracers));
 
  if (NumFiles < OMPcores)
     NC = NumFiles;
  else
     NC = OMPcores;  

  #pragma omp parallel for default(none) schedule(static) num_threads(NC) \
   private(i,snapshot,f1,f2,f3,Np,Header,SkipBlock,pos,vel,id,j,k,dummy)  \
   shared(Tracer,ScalePos,ScaleVel,stdout,NumFiles,FileTracers)  
  
  for (i=0; i<NumFiles; i++) {

      if (NumFiles == 1) 
         sprintf(snapshot,"%s",FileTracers);	  
      else 
         sprintf(snapshot,"%s.%d",FileTracers,i);	  

      f1 = SafeOpen(snapshot,"r"); // Pos     
      f2 = SafeOpen(snapshot,"r"); // Vel     
      f3 = SafeOpen(snapshot,"r"); // ID

      fread(&dummy,sizeof(int),1,f1);
      fread(&Header,sizeof(struct GadgetHeader),1,f1);
      fread(&dummy,sizeof(int),1,f1);
      fread(&dummy,sizeof(int),1,f1);

      Np = Header.Npart[1];

      SkipBlock = sizeof(int) + sizeof(struct GadgetHeader) + sizeof(int) 
                + sizeof(int) + 3*sizeof(float)*Np + 2*sizeof(int);
      fseek(f2,SkipBlock,SEEK_CUR);     

      SkipBlock = sizeof(int) + sizeof(struct GadgetHeader) + sizeof(int) 
                + sizeof(int) + 3*sizeof(float)*Np + sizeof(int)
                + sizeof(int) + 3*sizeof(float)*Np + 2*sizeof(int);
      fseek(f3,SkipBlock,SEEK_CUR);   

      for (j=0; j<Np; j++) {

	  fread(pos,sizeof(float),3,f1);    
	  fread(vel,sizeof(float),3,f2);    
	  fread(&id,sizeof(float),1,f3);  
	  id--;

	  for (k=0; k<3; k++) {
	      Tracer[id].Pos[k] = pos[k]*ScalePos; 
	      Tracer[id].Vel[k] = vel[k]*sqrt(Header.Time)*ScaleVel;
	  }

      }

      fclose(f1);
      fclose(f2);
      fclose(f3);
  }
  

}

void ReadTracers_GADGET2()
{

  struct GadgetHeader {
    int      Npart[6];
    double   Mass[6];
    double   Time;
    double   Redshift;
    int      Flag_1;
    int      Flag_2;
    int      NpartTotal[6];
    int      Flag_3;
    int      NumFiles;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam; 
    char     fill[96];  
  } Header;

  int   i,j,k,SkipBlock,SizeBlock,dummy,Np,NC;
  int   *Nstart,*NpFile,id;
  float pos[3],vel[3];
  FILE  *f1,*f2;
  char  snapshot[MAXCHAR],key[5],buffer[5];

  Nstart = (int *) malloc(NumFiles*sizeof(int));
  NpFile = (int *) malloc(NumFiles*sizeof(int));

  for (i=0; i<NumFiles; i++) {

      if (NumFiles == 1) 
         sprintf(snapshot,"%s",FileTracers);	  
      else 
         sprintf(snapshot,"%s.%d",FileTracers,i);	  

      f1 = SafeOpen(snapshot,"r");

      strcpy(key,"HEAD\0");
      do {
         fread(&dummy,sizeof(int),1,f1);
         fread(buffer,4*sizeof(char),1,f1); buffer[4]='\0';
         fread(&SizeBlock,sizeof(int),1,f1);
         if (strcmp(key,buffer) == 0) break;
         fseek(f1,SizeBlock+sizeof(int),SEEK_CUR);	  
      } while(1);
      fread(&dummy,sizeof(int),1,f1);
      fread(&dummy,sizeof(int),1,f1);
      fread(&Header,sizeof(struct GadgetHeader),1,f1); 
      fclose(f1);

      NpFile[i] = Header.Npart[1];
  }

  NumTrac = Header.NpartTotal[1];
  
  if (Header.BoxSize*ScalePos != BoxSize || Header.NumFiles != NumFiles) {
     fprintf(stdout,"\nError. Missmatch with Gadget header.\n");
     fprintf(stdout,"BoxSize = %f (%f in inputfile)\n",Header.BoxSize*ScalePos,BoxSize);
     fprintf(stdout,"NumFiles = %d (%d in inputfile)\n",Header.NumFiles,NumFiles);
     fflush(stdout);
     exit(EXIT_FAILURE);
  }

  if ((RSDist == 1 || GDist == 1) && (Header.Omega0 != OmegaMatter || 
       Header.OmegaLambda != OmegaLambda || Header.Redshift != Redshift || 
       Header.HubbleParam*100.0 != Hubble)) { 
     fprintf(stdout,"\nError. Missmatch with Gadget header.\n");
     fprintf(stdout,"OmegaMatter = %f (%f in inputfile)\n",Header.Omega0,OmegaMatter);
     fprintf(stdout,"OmegaLambda = %f (%f in inputfile)\n",Header.OmegaLambda,OmegaLambda);
     fprintf(stdout,"Hubble = %f (%f in inputfile)\n",Header.HubbleParam*100.0,Hubble);
     fprintf(stdout,"Redshift = %f (%f in inputfile)\n",Header.Redshift,Redshift);
     fflush(stdout);
     exit(EXIT_FAILURE);
  }
  
  Tracer = (struct tracers *) malloc(NumTrac*sizeof(struct tracers));
 
  for (i=0; i<NumFiles; i++) {
      Nstart[i] = 0;
      for (j=0; j<i; j++)
	  Nstart[i] += NpFile[j];      
  }


  if (NumFiles < OMPcores)
     NC = NumFiles;
  else
     NC = OMPcores;  

  #pragma omp parallel for default(none) schedule(static) num_threads(NC) \
   private(i,snapshot,f1,f2,Np,Header,SkipBlock,pos,vel,id,j,k,dummy,key, \
           buffer,SizeBlock)  \
   shared(Tracer,ScalePos,ScaleVel,stdout,NumFiles,FileTracers,Nstart)  
  
  for (i=0; i<NumFiles; i++) {

      if (NumFiles == 1) 
         sprintf(snapshot,"%s",FileTracers);	  
      else 
         sprintf(snapshot,"%s.%d",FileTracers,i);	  

      f1 = SafeOpen(snapshot,"r"); // Pos     
      f2 = SafeOpen(snapshot,"r"); // Vel     

      strcpy(key,"HEAD\0");
      do {
         fread(&dummy,sizeof(int),1,f1);
         fread(buffer,4*sizeof(char),1,f1); buffer[4]='\0';
         fread(&SizeBlock,sizeof(int),1,f1);
         if (strcmp(key,buffer) == 0) break;
         fseek(f1,SizeBlock+sizeof(int),SEEK_CUR);	  
      } while(1);
      fread(&dummy,sizeof(int),1,f1);
      fread(&dummy,sizeof(int),1,f1);
      fread(&Header,sizeof(struct GadgetHeader),1,f1);
      fread(&dummy,sizeof(int),1,f1);

      Np = Header.Npart[1];

      strcpy(key,"POS \0");
      do {
         fread(&dummy,sizeof(int),1,f1);
         fread(buffer,4*sizeof(char),1,f1); buffer[4]='\0';
         fread(&SizeBlock,sizeof(int),1,f1);
         if (strcmp(key,buffer) == 0) break;
         fseek(f1,SizeBlock+sizeof(int),SEEK_CUR);	  
      } while(1);
      fread(&dummy,sizeof(int),1,f1);
      fread(&dummy,sizeof(int),1,f1);

      strcpy(key,"VEL \0");
      do {
         fread(&dummy,sizeof(int),1,f2);
         fread(buffer,4*sizeof(char),1,f2); buffer[4]='\0';
         fread(&SizeBlock,sizeof(int),1,f2);
         if (strcmp(key,buffer) == 0) break;
         fseek(f2,SizeBlock+sizeof(int),SEEK_CUR);	  
      } while(1);
      fread(&dummy,sizeof(int),1,f2);
      fread(&dummy,sizeof(int),1,f2);

      //fseek(f1,sizeof(float)*Header.Npart[1]*3,SEEK_CUR);    
      //fseek(f2,sizeof(float)*Header.Npart[1]*3,SEEK_CUR);

      id = Nstart[i];
      for (j=0; j<Np; j++) {

	  fread(pos,sizeof(float),3,f1);    
	  fread(vel,sizeof(float),3,f2);    

	  for (k=0; k<3; k++) {
	      Tracer[id].Pos[k] = pos[k]*ScalePos; 
	      Tracer[id].Vel[k] = vel[k]*sqrt(Header.Time)*ScaleVel;
	  }

	  id++;
      }

      fclose(f1);
      fclose(f2);
  }

  free(Nstart);
  free(NpFile);

}

void ReadTracers_MXXL()
{

   int     i,j,k,NumTot,N,id;
   float   Pos[3],Vel[3];
   char    filename[MAXCHAR],basename[MAXCHAR];
   FILE    *fd;
	
   sprintf(basename,"%s/z%4.2f/halos_z%4.2f_part",FileTracers,Redshift,Redshift);
   fprintf(logfile," | Files = %s (%d files) \n",basename,NumFiles);

   for (i=1; i<=NumFiles; i++) {
       sprintf(filename,"%s%02d",basename,i);
       fd = SafeOpen(filename,"r");
       fread(&N,sizeof(int),1,fd);
       fread(&NumTot,sizeof(int),1,fd);
       NumTrac += N;
       fclose(fd);
   }

   if (NumTrac != NumTot) {
      fprintf(stdout,"Error. Missing halos: NumTrac = %d, NumTot = %d\n",NumTrac,NumTot);
      fflush(stdout);
      exit(EXIT_FAILURE);
   }

   Tracer = (struct tracers *) malloc(NumTrac*sizeof(struct tracers));

   for (j=1; j<=NumFiles; j++) { 
       sprintf(filename,"%s%02d",basename,j);
       fd = SafeOpen(filename,"r");
       fread(&N,sizeof(int),1,fd);
       fread(&NumTot,sizeof(int),1,fd);
       for (i=0; i<N; i++) {
	   fread(&id,sizeof(int),1,fd);
	   fread(Pos,sizeof(float),3,fd);
	   fread(Vel,sizeof(float),3,fd);

	   for (k=0; k<3; k++) {
	       Tracer[id].Pos[k] = Pos[k]*ScalePos;
	       Tracer[id].Vel[k] = Vel[k]*ScaleVel;
	   }
       }
       fclose(fd);
   }

}

void RedshiftSpaceDistortions()
{

   int    i;
   double Hubble_z;
   double RSDFactor;	
   struct cosmology C = {OmegaMatter,(1.0 - OmegaMatter - OmegaLambda),OmegaLambda};

   fprintf(logfile," | Applying redshift-space distortions: LOS = z-axis, POS = xy-plane\n");

   if (Redshift == 0.0) {
     RSDFactor = 1.0/100.0;	   
   } else { 
     Hubble_z = 100.0*E(Redshift,&C);
     RSDFactor = (1.0 + Redshift)/Hubble_z;
   }

   fprintf(logfile," | RSD factor = (1+z)/H(z) = %f [h⁻¹Mpc/(km/s)]\n",RSDFactor);
   fflush(logfile);

   for (i=0; i<NumTrac; i++) { 
       Tracer[i].Pos[2] += Tracer[i].Vel[2]*RSDFactor;
       if (Tracer[i].Pos[2] < 0.0    ) Tracer[i].Pos[2] += LBox[2];
       if (Tracer[i].Pos[2] > LBox[2]) Tracer[i].Pos[2] -= LBox[2];        
   }	       

}

void GeometricalDistortions()
{

   int    i;
   double Hubble_z,Distance_z;	
   double FidHubble_z,FidDistance_z;	
   double GDFactor_LOS,GDFactor_POS;
   
   struct cosmology C = {OmegaMatter,
	                 OmegaLambda,
			 1.0 - OmegaMatter - OmegaLambda,
                         Hubble};
   struct cosmology FC = {FidOmegaMatter,
	                  FidOmegaLambda,
			  1.0 - FidOmegaMatter - FidOmegaLambda,
                          FidHubble};

   fprintf(logfile," | Applying fiducial cosmology distortions: LOS = z-axis, POS = xy-plane\n");
   fprintf(logfile," | True cosmology: (OmM,OmL,OmK,H0) = (%4.2f,%4.2f,%4.2f,%4.2f)\n",C.OmM,C.OmL,C.OmK,C.Hub);

   Hubble_z = C.Hub*E(Redshift,&C);
   Distance_z = AngularDistance(Redshift,&C);
   
   fprintf(logfile," | Fiducial cosmology: (OmM,OmL,OmK,H0) = (%4.2f,%4.2f,%4.2f,%4.2f)\n",FC.OmM,FC.OmL,FC.OmK,FC.Hub);

   FidHubble_z = FC.Hub*E(Redshift,&FC);
   FidDistance_z = AngularDistance(Redshift,&FC);

   GDFactor_LOS = Hubble_z/FidHubble_z;
   GDFactor_POS = FidDistance_z/Distance_z;

   fprintf(logfile," | GD factor: LOS = H(z)/FidH(z) = %f \n",GDFactor_LOS);
   fprintf(logfile," |            POS = FidD(z)/D(z) = %f \n",GDFactor_POS);
   fflush(logfile);

   LBox[0] *= GDFactor_POS;
   LBox[1] *= GDFactor_POS;
   LBox[2] *= GDFactor_LOS;

   for (i=0; i<NumTrac; i++) { 
       Tracer[i].Pos[0] *= GDFactor_POS;
       Tracer[i].Pos[1] *= GDFactor_POS;
       Tracer[i].Pos[2] *= GDFactor_LOS;
   }	

}

void WriteVoids()
{
   int     i;
   FILE    *fd;
   clock_t t;

   fprintf(logfile,"\n WRITTING VOID CATALOGUE \n");
   t = clock();
   
   fd = SafeOpen(FileVoids,"w");

   float dx[3],disp;
   int   k;

   for (i=0; i<NumVoid; i++) {
       if (Void[i].ToF) {

          fprintf(fd," %8.5f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %5d\n",
	  	       Void[i].Rad,Void[i].Pos[0],Void[i].Pos[1],Void[i].Pos[2],Void[i].Vel[0],Void[i].Vel[1],
	  	       Void[i].Vel[2],Void[i].Delta,Void[i].Dtype,Void[i].Poisson,Void[i].Nran);   
       }
   }
   fclose(fd);

   StepName.push_back("Writting void catalogue");
   StepTime.push_back(Time(t,1));

}

