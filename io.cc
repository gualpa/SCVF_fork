
#include "allvars.h"
#include "io.h"
#include "tools.h"
#include "cosmology.h"

void read_input_file(char *filename)
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
  addr[nt] = &VarConfig.BoxSize;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"MaxRadiusSearch");
  addr[nt] = &VarConfig.MaxRadiusSearch;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"ProxyGridSize");
  addr[nt] = &VarConfig.ProxyGridSize;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"DeltaThreshold");
  addr[nt] = &VarConfig.DeltaThreshold;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"DeltaSeed");
  addr[nt] = &VarConfig.DeltaSeed;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"OverlapTol");
  addr[nt] = &VarConfig.OverlapTol;
  id[nt++] = DOUBLE;  

  strcpy(tag[nt],"FormatTracers");
  addr[nt] = &VarConfig.FormatTracers;
  id[nt++] = INT;

  strcpy(tag[nt],"NumFiles");
  addr[nt] = &VarConfig.NumFiles;
  id[nt++] = INT;

  strcpy(tag[nt],"FileTracers");
  addr[nt] = VarConfig.FileTracers;
  id[nt++] = STRING;

  strcpy(tag[nt],"FileVoids");
  addr[nt] = VarConfig.FileVoids;
  id[nt++] = STRING;

  strcpy(tag[nt],"ScalePos");
  addr[nt] = &Scale.Pos;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt],"ScaleVel");
  addr[nt] = &Scale.Vel;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"OMPcores");
  addr[nt] = &VarConfig.OMPcores;
  id[nt++] = INT;

  strcpy(tag[nt],"RadIncrement");
  addr[nt] = &VarConfig.RadIncrement;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt],"NumRanWalk");
  addr[nt] = &VarConfig.NumRanWalk;
  id[nt++] = INT;

  strcpy(tag[nt],"FracRadius");
  addr[nt] = &VarConfig.FracRadius;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"RSDist");
  addr[nt] = &VarConfig.RSDist;
  id[nt++] = INT;

  strcpy(tag[nt],"Redshift");
  addr[nt] = &VarConfig.Redshift;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"OmegaMatter");
  addr[nt] = &VarConfig.OmegaMatter;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"OmegaLambda");
  addr[nt] = &VarConfig.OmegaLambda;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"Hubble");
  addr[nt] = &VarConfig.Hubble;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"GDist");
  addr[nt] = &VarConfig.GDist;
  id[nt++] = INT;

  strcpy(tag[nt],"FidOmegaMatter");
  addr[nt] = &VarConfig,FidOmegaMatter;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"FidOmegaLambda");
  addr[nt] = &VarConfig.FidOmegaLambda;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"FidHubble");
  addr[nt] = &VarConfig.FidHubble;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"WriteProfiles");
  addr[nt] = &VarConfig.WriteProfiles;
  id[nt++] = INT;

  strcpy(tag[nt],"MinProfileDist");
  addr[nt] = &VarConfig.MinProfileDist;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"MaxProfileDist");
  addr[nt] = &VarConfig.MaxProfileDist;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"NumProfileBins");
  addr[nt] = &VarConfig.NumProfileBins;
  id[nt++] = INT;

  strcpy(tag[nt],"PathProfiles");
  addr[nt] = VarConfig.PathProfiles;
  id[nt++] = STRING;

  strcpy(tag[nt],"InnerShellVel");
  addr[nt] = &VarConfig.InnerShellVel;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"OuterShellVel");
  addr[nt] = &VarConfig.OuterShellVel;
  id[nt++] = DOUBLE; 

  fd = safe_open(filename,"r");

  if (VarConfig.RunFlag == 0) 
     sprintf(fname,"%s.log",filename);	  
  else
     sprintf(fname,"%s_%d.log",filename,VarConfig.RunFlag);	  

  VarConfig.logfile = safe_open(fname,"w");
  fprintf(VarConfig.logfile,"\n CONFIGURATION PARAMETERS USED \n\n");

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
    	    fprintf(VarConfig.logfile, " %-35s%g\n", buf1, *((double *) addr[j]));
    	    break;

    	    case STRING:
    	    strcpy((char *)addr[j], buf2);
    	    fprintf(VarConfig.logfile, " %-35s%s\n", buf1, buf2);
    	    break;

    	    case INT:
    	    *((int *) addr[j]) = atoi(buf2);
    	    fprintf(VarConfig.logfile, " %-35s%d\n", buf1, *((int *) addr[j]));
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

  fflush(VarConfig.logfile);

#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
}

void read_tracers()
{

   fprintf(VarConfig.logfile,"\n READING TRACERS \n");
   clock_t t = clock();

   switch (VarConfig.FormatTracers) {

      case 0: 
      fprintf(VarConfig.logfile," | Reading ASCII format \n");	   
      read_tracers_ascii();
      break;

      case 1:
      fprintf(VarConfig.logfile," | Reading GADGET \n");	   
      read_tracers_gadget();
      break;

      case 2:
      fprintf(VarConfig.logfile," | Reading MXXL format \n");	   
      read_tracers_mxxl();
      break;

      case 3:
      fprintf(VarConfig.logfile," | Reading BINARY format \n");	   
      read_tracers_binary();

   }

   float xmax = 0.0;
   float ymax = 0.0;
   float zmax = 0.0;
   float diff = 0.999;
   for (int i=0; i<VarConfig.NumTrac; i++) {
       if (Tracer[i].Pos[0] > xmax) xmax = Tracer[i].Pos[0];	   
       if (Tracer[i].Pos[1] > ymax) ymax = Tracer[i].Pos[1];	   
       if (Tracer[i].Pos[2] > zmax) zmax = Tracer[i].Pos[2];	   
   }
   if (xmax/BoxSize < diff || ymax/VarConfig.BoxSize < diff || zmax/VarConfig.BoxSize < diff) {
      fprintf(stdout,"\n Error. Wrong BoxSize? - MAX = (%f,%f,%f) \n",xmax,ymax,zmax);
      fflush(stdout);
      exit(EXIT_FAILURE);	
   } 

   VarConfig.LBox[0] = VarConfig.BoxSize;
   VarConfig.LBox[1] = VarConfig.BoxSize;
   VarConfig.LBox[2] = VarConfig.BoxSize;

   if (VarConfig.RSDist == 1) redshift_space_distortions();
   if (VarConfig.GDist == 1) geometrical_distortions();

   double Volume = VarConfig.LBox[0]*VarConfig.LBox[1]*VarConfig.LBox[2];
   VarConfig.MeanNumTrac = (double)VarConfig.NumTrac/Volume;
   VarConfig.MeanSeparation = cbrt(Volume/(double)VarConfig.NumTrac);
 
   fprintf(VarConfig.logfile," | Number of tracers = %d \n",VarConfig.NumTrac);
   fprintf(VarConfig.logfile," | Size of the box: x-axis = %f \n",VarConfig.LBox[0]);
   fprintf(VarConfig.logfile," |                  y-axis = %f \n",VarConfig.LBox[1]);
   fprintf(VarConfig.logfile," |                  z-axis = %f \n",VarConfig.LBox[2]);
   fprintf(VarConfig.logfile," | Mean number density [h³/Mpc³] = %e \n",VarConfig.MeanNumTrac);
   fprintf(VarConfig.logfile," | Mean separation [Mpc/h] = %e \n",VarConfig.MeanSeparation);

   VarConfig.StepName.push_back("Reading tracers");
   VarConfig.StepTime.push_back(get_time(t,1));

}

void read_tracers_ascii()
{
   
   VarConfig.NumTrac = 0;	
   int NumTot = count_lines(VarConfig.FileTracers);
   FILE *fd = safe_open(VarConfig.FileTracers,"r");
   float dummy ; // agrego Seba par poder levantar archivo de 7 columnas
   for (int i=0; i<NumTot; i++) {

       //if (NumTrac > 250000) break	   

       Tracer.push_back(tracers());

       fscanf(fd,"%f %f %f %f %f %f %f\n",&Tracer[i].Pos[0],&Tracer[i].Pos[1],&Tracer[i].Pos[2],	   
                                        &Tracer[i].Vel[0],&Tracer[i].Vel[1],&Tracer[i].Vel[2],&dummy);
       Tracer[VarConfig.NumTrac].Pos[0] *= Scale.Pos;
       Tracer[VarConfig.NumTrac].Pos[1] *= Scale.Pos;
       Tracer[VarConfig.NumTrac].Pos[2] *= Scale.Pos;
       Tracer[VarConfig.NumTrac].Vel[0] *= Scale.Vel;
       Tracer[VarConfig.NumTrac].Vel[1] *= Scale.Vel;
       Tracer[VarConfig.NumTrac].Vel[2] *= Scale.Vel;
       VarConfig.NumTrac++;
   }
   fclose(fd);

}

void read_tracers_binary()
{
   FILE *fd = safe_open(VarConfig.FileTracers,"r");
   int NumTot;

   VarConfig.NumTrac = 0;
   fread(&NumTot,sizeof(int),1,fd);

   for (int i=0; i<NumTot; i++) {
       Tracer.push_back(tracers());

       fread(&Tracer[VarConfig.NumTrac].Pos[0],sizeof(float),3,fd);
       fread(&Tracer[VarConfig.NumTrac].Vel[0],sizeof(float),3,fd);
       
       Tracer[VarConfig.NumTrac].Pos[0] *= Scale.Pos;
       Tracer[VarConfig.NumTrac].Pos[1] *= Scale.Pos;
       Tracer[VarConfig.NumTrac].Pos[2] *= Scale.Pos;
       Tracer[VarConfig.NumTrac].Vel[0] *= Scale.Vel;
       Tracer[VarConfig.NumTrac].Vel[1] *= Scale.Vel;
       Tracer[VarConfig.NumTrac].Vel[2] *= Scale.Vel;
       VarConfig.NumTrac++;
   }
   fclose(fd);

}

void read_tracers_gadget()
{
#define SKIP fread(&dummy,sizeof(int),1,f1)

/*
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
*/

  struct GadgetHeader {
    long long Npart[2];
    long long NpartTotal[2];
    double    Mass[2];
    double    Time;
    double    Redshift;
    double    BoxSize;
    int       NumFiles;
    long long Ntrees;
    long long NtreesTot;
  } Header;

  typedef float postype;
  typedef int idtype;

  int  i,j,k,dummy,Np,NC,jump;
  int  possize,velsize,idsize;
  FILE *f1;
  idtype *id;
  postype *pos,*vel;
  char snapshot[MAXCHAR];

  if (VarConfig.NumFiles == 1) 
     sprintf(snapshot,"%s",VarConfig.FileTracers);	  
  else 
     sprintf(snapshot,"%s.0",VarConfig.FileTracers);	  

  f1 = safe_open(snapshot,"r");

  SKIP; fread(&Header,sizeof(struct GadgetHeader),1,f1); SKIP; 

  fread(&possize,sizeof(int),1,f1);
  fseek(f1,possize+sizeof(int),SEEK_CUR);
  fread(&velsize,sizeof(int),1,f1);
  fseek(f1,velsize+sizeof(int),SEEK_CUR);
  fread(&idsize,sizeof(int),1,f1);
  
  possize /= (3*Header.Npart[1]); 
  velsize /= (3*Header.Npart[1]); 
  idsize /= Header.Npart[1]; 

  fclose(f1);

  if (Header.BoxSize*Scale.Pos != VarConfig.BoxSize || Header.NumFiles != VarConfig.NumFiles) {
     fprintf(stdout,"\nError. Missmatch with Gadget header.\n");
     fprintf(stdout,"BoxSize = %f (%f in inputfile)\n",Header.BoxSize*Scale.Pos,VarConfig.BoxSize);
     fprintf(stdout,"NumFiles = %d (%d in inputfile)\n",Header.NumFiles,VarConfig.NumFiles);
     fflush(stdout);
     exit(EXIT_FAILURE);
  }

  if ((VarConfig.RSDist == 1 || VarConfig.GDist == 1) && Header.Redshift != VarConfig.Redshift) { 
     fprintf(stdout,"\nError. Missmatch with Gadget header.\n");
     fprintf(stdout,"Redshift = %f (%f in inputfile)\n",Header.Redshift,Redshift);
     fflush(stdout);
     exit(EXIT_FAILURE);
  }

  if (possize != sizeof(postype) || velsize != sizeof(postype) || idsize != sizeof(idtype)) {
     fprintf(stdout,"\nError. Missmatch data types in Gadget file.\n");
     fprintf(stdout,"Pos: defined as %d bytes, in file %d bytes. \n",sizeof(postype),possize);
     fprintf(stdout,"Vel: defined as %d bytes, in file %d bytes. \n",sizeof(postype),velsize);
     fprintf(stdout,"IDs: defined as %d bytes, in file %d bytes. \n",sizeof(idtype),idsize);
     fflush(stdout);
     exit(EXIT_FAILURE);
  } 

  VarConfig.NumTrac = Header.NpartTotal[1];
  for (int i=0; i<VarConfig.NumTrac; i++) Tracer.push_back(tracers());

  if (VarConfig.NumFiles < VarConfig.OMPcores)
     NC = VarConfig.NumFiles;
  else
     NC = VarConfig.OMPcores;  

  #pragma omp parallel for default(none) schedule(static) num_threads(NC) \
   private(i,snapshot,f1,Np,Header,pos,vel,id,j,k,dummy)  \
   shared(Tracer,Scale.Pos,Scale.Vel,stdout,VarConfig.NumFiles,FileTracers)  
  
  for (i=0; i<VarConfig.NumFiles; i++) {

      if (VarConfig.NumFiles == 1) 
         sprintf(snapshot,"%s",FileTracers);	  
      else 
         sprintf(snapshot,"%s.%d",FileTracers,i);	  

      f1 = safe_open(snapshot,"r");     

      SKIP; fread(&Header,sizeof(struct GadgetHeader),1,f1); SKIP;  

      Np = Header.Npart[1];
      pos = (postype *) malloc(3*Np*sizeof(postype));
      vel = (postype *) malloc(3*Np*sizeof(postype));
      id = (idtype *) malloc(Np*sizeof(idtype));

      SKIP; fread(pos,sizeof(postype),3*Np,f1); SKIP;
      SKIP; fread(vel,sizeof(postype),3*Np,f1); SKIP;
      SKIP; fread(id,sizeof(idtype),Np,f1); SKIP;

      fclose(f1);

      for (j=0; j<Np; j++) {
	  for (k=0; k<3; k++) {
	      Tracer[id[j]-1].Pos[k] = pos[3*j+k]*Scale.Pos; 
	      Tracer[id[j]-1].Vel[k] = vel[3*j+k]*sqrt(Header.Time)*Scale.Vel;
	  }
      }

      free(pos);
      free(vel);
      free(id);
  }
#undef SKIP 
}

void read_tracers_mxxl()
{

   int     i,j,k,NumTot,N,id;
   float   Pos[3],Vel[3];
   char    filename[MAXCHAR],basename[MAXCHAR];
   FILE    *fd;
	
   sprintf(basename,"%s/z%4.2f/halos_z%4.2f_part",VarConfig.FileTracers,VarConfig.Redshift,VarConfig.Redshift);
   fprintf(VarConfig.logfile," | Files = %s (%d files) \n",basename,VarConfig.NumFiles);

   for (i=1; i<=VarConfig.NumFiles; i++) {
       sprintf(filename,"%s%02d",basename,i);
       fd = safe_open(filename,"r");
       fread(&N,sizeof(int),1,fd);
       fread(&NumTot,sizeof(int),1,fd);
       VarConfig.NumTrac += N;
       fclose(fd);
   }

   if (VarConfig.NumTrac != NumTot) {
      fprintf(stdout,"Error. Missing halos: NumTrac = %d, NumTot = %d\n",VarConfig.NumTrac,NumTot);
      fflush(stdout);
      exit(EXIT_FAILURE);
   }

   for (int i=0; i<VarConfig.NumTrac; i++) Tracer.push_back(tracers());

   for (j=1; j<=VarConfig.NumFiles; j++) { 
       sprintf(filename,"%s%02d",basename,j);
       fd = safe_open(filename,"r");
       fread(&N,sizeof(int),1,fd);
       fread(&NumTot,sizeof(int),1,fd);
       for (i=0; i<N; i++) {
	   fread(&id,sizeof(int),1,fd);
	   fread(Pos,sizeof(float),3,fd);
	   fread(Vel,sizeof(float),3,fd);

	   for (k=0; k<3; k++) {
	       Tracer[id].Pos[k] = Pos[k]*Scale.Pos;
	       Tracer[id].Vel[k] = Vel[k]*Scale.Vel;
	   }
       }
       fclose(fd);
   }

}

void redshift_space_distortions()
{

   int    i;
   double Hubble_z;
   double RSDFactor;	
   struct cosmoparam C = {
      VarConfig.OmegaMatter,VarConfig.OmegaLambda,
      (1.0 - VarConfig.OmegaMatter - VarConfig.OmegaLambda),
      VarConfig.Hubble
   };

   fprintf(VarConfig.logfile," | Applying redshift-space distortions: LOS = z-axis, POS = xy-plane\n");

   if (VarConfig.Redshift == 0.0) {
     RSDFactor = 1.0/100.0;	   
   } else { 
     Hubble_z = 100.0*evolution_param(VarConfig.Redshift,&C);
     RSDFactor = (1.0 + VarConfig.Redshift)/Hubble_z;
   }

   fprintf(VarConfig.logfile," | RSD factor = (1+z)/H(z) = %f [h⁻¹Mpc/(km/s)]\n",RSDFactor);
   fflush(VarConfig.logfile);

   for (i=0; i<VarConfig.NumTrac; i++) { 
       Tracer[i].Pos[2] += Tracer[i].Vel[2]*RSDFactor;
       if (Tracer[i].Pos[2] < 0.0    ) Tracer[i].Pos[2] += VarConfig.LBox[2];
       if (Tracer[i].Pos[2] > VarConfig.LBox[2]) Tracer[i].Pos[2] -= VarConfig.LBox[2];        
   }	       

}

void geometrical_distortions()
{

   int    i;
   double Hubble_z,Distance_z;	
   double FidHubble_z,FidDistance_z;	
   double GDFactor_LOS,GDFactor_POS;
   struct cosmoparam C = {
      VarConfig.OmegaMatter,
      VarConfig.OmegaLambda,
		1.0 - VarConfig.OmegaMatter - VarConfig.OmegaLambda,
      VarConfig.Hubble
   };
   struct cosmoparam FC = {
      VarConfig.FidOmegaMatter,
	   VarConfig.FidOmegaLambda,
		1.0 - VarConfig.FidOmegaMatter - VarConfig.FidOmegaLambda,
      VarConfig.FidHubble
   };

   fprintf(VarConfig.logfile," | Applying fiducial cosmology distortions: LOS = z-axis, POS = xy-plane\n");
   fprintf(VarConfig.logfile," | True cosmology: (OmM,OmL,OmK,H0) = (%4.2f,%4.2f,%4.2f,%4.2f)\n",C.OmM,C.OmL,C.OmK,C.Hub);

   Hubble_z = C.Hub*evolution_param(VarConfig.Redshift,&C);
   Distance_z = angular_distance(VarConfig.Redshift,&C);
   
   fprintf(VarConfig.logfile," | Fiducial cosmology: (OmM,OmL,OmK,H0) = (%4.2f,%4.2f,%4.2f,%4.2f)\n",FC.OmM,FC.OmL,FC.OmK,FC.Hub);

   FidHubble_z = FC.Hub*evolution_param(Redshift,&FC);
   FidDistance_z = angular_distance(Redshift,&FC);

   GDFactor_LOS = Hubble_z/FidHubble_z;
   GDFactor_POS = FidDistance_z/Distance_z;

   fprintf(VarConfig.logfile," | GD factor: LOS = H(z)/FidH(z) = %f \n",GDFactor_LOS);
   fprintf(VarConfig.logfile," |            POS = FidD(z)/D(z) = %f \n",GDFactor_POS);
   fflush(VarConfig.logfile);

   VarConfig.LBox[0] *= GDFactor_POS;
   VarConfig.LBox[1] *= GDFactor_POS;
   VarConfig.LBox[2] *= GDFactor_LOS;

   for (i=0; i<VarConfig.NumTrac; i++) { 
       Tracer[i].Pos[0] *= GDFactor_POS;
       Tracer[i].Pos[1] *= GDFactor_POS;
       Tracer[i].Pos[2] *= GDFactor_LOS;
   }	

}

void write_voids()
{
   int     i;
   FILE    *fd;
   clock_t t;

   fprintf(VarConfig.logfile,"\n WRITTING VOID CATALOGUE \n");
   t = clock();
   
   fd = safe_open(VarConfigFileVoids,"w");

   for (i=0; i<VarConfig.NumVoid; i++) {
       if (Void[i].ToF) {

          fprintf(fd," %8.5f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %5d\n",
	  	       Void[i].Rad,Void[i].Pos[0],Void[i].Pos[1],Void[i].Pos[2],Void[i].Vel[0],Void[i].Vel[1],
	  	       Void[i].Vel[2],Void[i].Delta,Void[i].Dtype,Void[i].Poisson,Void[i].Nran);   
       }
   }
   fclose(fd);

   VarConfig.StepName.push_back("Writting void catalogue");
   VarConfig.StepTime.push_back(get_time(t,1));

}

/*
void read_tracers_gadget2_format2()
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

      f1 = safe_open(snapshot,"r");

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

  if (Header.BoxSize*Scale.Pos != BoxSize || Header.NumFiles != NumFiles) {
     fprintf(stdout,"\nError. Missmatch with Gadget header.\n");
     fprintf(stdout,"BoxSize = %f (%f in inputfile)\n",Header.BoxSize*Scale.Pos,BoxSize);
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
  
  NumTrac = Header.NpartTotal[1];
  for (int i=0; i<NumTrac; i++) Tracer.push_back(tracers());
 
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
   shared(Tracer,Scale.Pos,Scale.Vel,stdout,NumFiles,FileTracers,Nstart)  
  
  for (i=0; i<NumFiles; i++) {

      if (NumFiles == 1) 
         sprintf(snapshot,"%s",FileTracers);	  
      else 
         sprintf(snapshot,"%s.%d",FileTracers,i);	  

      f1 = safe_open(snapshot,"r"); // Pos     
      f2 = safe_open(snapshot,"r"); // Vel     

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
	      Tracer[id].Vel[k] = vel[k]*sqrt(Header.Time)*Scale.Vel;
	  }

	  id++;
      }

      fclose(f1);
      fclose(f2);
  }

  free(Nstart);
  free(NpFile);

}
 */
