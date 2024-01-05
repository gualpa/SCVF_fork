
#include "allvars.h"
#include "io.h"
#include "tools.h"
#include "cosmology.h"

/**
 * @brief Carga los valores del archivo de configuración.
 *
 * Esta función lee el archivo de configuración y lo carga en la variable que devuelve.
 *
 * @param filename nombre del archivo de configuración
 * @param VarConfigAux estructura donde estan almacenados valores de configuración del programa.
 * @return Una estructura con los valores de configuración.
 */
varConfiguration read_input_file(char *filename, varConfiguration VarConfigAux, logs &LogAux)
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
  scale_struct ScaleAux;
  nt = 0;

  strcpy(tag[nt],"BoxSize");
  addr[nt] = &VarConfigAux.BoxSize;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"MaxRadiusSearch");
  addr[nt] = &VarConfigAux.MaxRadiusSearch;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"ProxyGridSize");
  addr[nt] = &VarConfigAux.ProxyGridSize;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"DeltaThreshold");
  addr[nt] = &VarConfigAux.DeltaThreshold;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"DeltaSeed");
  addr[nt] = &VarConfigAux.DeltaSeed;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"OverlapTol");
  addr[nt] = &VarConfigAux.OverlapTol;
  id[nt++] = DOUBLE;  

  strcpy(tag[nt],"FormatTracers");
  addr[nt] = &VarConfigAux.FormatTracers;
  id[nt++] = INT;

  strcpy(tag[nt],"NumFiles");
  addr[nt] = &VarConfigAux.NumFiles;
  id[nt++] = INT;

  strcpy(tag[nt],"FileTracers");
  addr[nt] = VarConfigAux.FileTracers;
  id[nt++] = STRING;

  strcpy(tag[nt],"FileVoids");
  addr[nt] = VarConfigAux.FileVoids;
  id[nt++] = STRING;

  strcpy(tag[nt],"ScalePos");
  addr[nt] = &ScaleAux.Pos;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt],"ScaleVel");
  addr[nt] = &ScaleAux.Vel;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"OMPcores");
  addr[nt] = &VarConfigAux.OMPcores;
  id[nt++] = INT;

  strcpy(tag[nt],"RadIncrement");
  addr[nt] = &VarConfigAux.RadIncrement;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt],"NumRanWalk");
  addr[nt] = &VarConfigAux.NumRanWalk;
  id[nt++] = INT;

  strcpy(tag[nt],"FracRadius");
  addr[nt] = &VarConfigAux.FracRadius;
  id[nt++] = DOUBLE;

  strcpy(tag[nt],"RSDist");
  addr[nt] = &VarConfigAux.RSDist;
  id[nt++] = INT;

  strcpy(tag[nt],"Redshift");
  addr[nt] = &VarConfigAux.Redshift;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"OmegaMatter");
  addr[nt] = &VarConfigAux.OmegaMatter;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"OmegaLambda");
  addr[nt] = &VarConfigAux.OmegaLambda;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"Hubble");
  addr[nt] = &VarConfigAux.Hubble;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"GDist");
  addr[nt] = &VarConfigAux.GDist;
  id[nt++] = INT;

  strcpy(tag[nt],"FidOmegaMatter");
  addr[nt] = &VarConfigAux.FidOmegaMatter;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"FidOmegaLambda");
  addr[nt] = &VarConfigAux.FidOmegaLambda;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"FidHubble");
  addr[nt] = &VarConfigAux.FidHubble;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"WriteProfiles");
  addr[nt] = &VarConfigAux.WriteProfiles;
  id[nt++] = INT;

  strcpy(tag[nt],"MinProfileDist");
  addr[nt] = &VarConfigAux.MinProfileDist;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"MaxProfileDist");
  addr[nt] = &VarConfigAux.MaxProfileDist;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"NumProfileBins");
  addr[nt] = &VarConfigAux.NumProfileBins;
  id[nt++] = INT;

  strcpy(tag[nt],"PathProfiles");
  addr[nt] = VarConfigAux.PathProfiles;
  id[nt++] = STRING;

  strcpy(tag[nt],"InnerShellVel");
  addr[nt] = &VarConfigAux.InnerShellVel;
  id[nt++] = DOUBLE; 

  strcpy(tag[nt],"OuterShellVel");
  addr[nt] = &VarConfigAux.OuterShellVel;
  id[nt++] = DOUBLE; 

  fd = safe_open(filename,"r");

  if (VarConfigAux.RunFlag == 0)
     sprintf(fname,"%s.log",filename);
  else
     sprintf(fname,"%s_%d.log",filename,VarConfigAux.RunFlag);

  LogAux.logfile = safe_open(fname,"w");
  fprintf(LogAux.logfile,"\n CONFIGURATION PARAMETERS USED \n\n");

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
          fprintf(LogAux.logfile, " %-35s%g\n", buf1, *((double *) addr[j]));
    	    break;

    	    case STRING:
    	    strcpy((char *)addr[j], buf2);
          fprintf(LogAux.logfile, " %-35s%s\n", buf1, buf2);
    	    break;

    	    case INT:
          *((int *) addr[j]) = atoi(buf2);
          fprintf(LogAux.logfile, " %-35s%d\n", buf1, *((int *) addr[j]));
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

  fflush(LogAux.logfile);

#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
  VarConfigAux.Scale = ScaleAux;
  return VarConfigAux;
}

void read_tracers(varConfiguration &VarConfigAux, logs &LogAux)
{
   clock_t t = clock();
   fprintf(LogAux.logfile,"\n READING TRACERS \n");

   switch (VarConfigAux.FormatTracers) {

      case 0: 
          fprintf(LogAux.logfile," | Reading ASCII format \n");
          read_tracers_ascii(VarConfigAux);
      break;

      case 1:
          fprintf(LogAux.logfile," | Reading GADGET \n");
          read_tracers_gadget(VarConfigAux);
      break;

      case 2:
          fprintf(LogAux.logfile," | Reading MXXL format \n");
          read_tracers_mxxl(VarConfigAux, LogAux);
      break;

      case 3:
          fprintf(LogAux.logfile," | Reading BINARY format \n");
          read_tracers_binary(VarConfigAux);

   }

   float xmax = 0.0;
   float ymax = 0.0;
   float zmax = 0.0;
   float diff = 0.999;
   for (int i=0; i<VarConfigAux.NumTrac; i++) {
       if (Tracer[i].Pos[0] > xmax) xmax = Tracer[i].Pos[0];	   
       if (Tracer[i].Pos[1] > ymax) ymax = Tracer[i].Pos[1];	   
       if (Tracer[i].Pos[2] > zmax) zmax = Tracer[i].Pos[2];	   
   }
   if (xmax/VarConfigAux.BoxSize < diff || ymax/VarConfigAux.BoxSize < diff || zmax/VarConfigAux.BoxSize < diff) {
      fprintf(stdout,"\n Error. Wrong BoxSize? - MAX = (%f,%f,%f) \n",xmax,ymax,zmax);
      fflush(stdout);
      exit(EXIT_FAILURE);	
   } 

   VarConfigAux.LBox[0] = VarConfigAux.BoxSize;
   VarConfigAux.LBox[1] = VarConfigAux.BoxSize;
   VarConfigAux.LBox[2] = VarConfigAux.BoxSize;

   if (VarConfigAux.RSDist == 1) redshift_space_distortions(VarConfigAux, LogAux);
   if (VarConfigAux.GDist == 1) geometrical_distortions(VarConfigAux, LogAux);

   double Volume = VarConfigAux.LBox[0]*VarConfigAux.LBox[1]*VarConfigAux.LBox[2];
   VarConfigAux.MeanNumTrac = (double)VarConfigAux.NumTrac/Volume;
   VarConfigAux.MeanSeparation = cbrt(Volume/(double)VarConfigAux.NumTrac);
 
   fprintf(LogAux.logfile," | Number of tracers = %d \n",VarConfigAux.NumTrac);
   fprintf(LogAux.logfile," | Size of the box: x-axis = %f \n",VarConfigAux.LBox[0]);
   fprintf(LogAux.logfile," |                  y-axis = %f \n",VarConfigAux.LBox[1]);
   fprintf(LogAux.logfile," |                  z-axis = %f \n",VarConfigAux.LBox[2]);
   fprintf(LogAux.logfile," | Mean number density [h³/Mpc³] = %e \n",VarConfigAux.MeanNumTrac);
   fprintf(LogAux.logfile," | Mean separation [Mpc/h] = %e \n",VarConfigAux.MeanSeparation);

   LogAux.StepName.push_back("Reading tracers");
   LogAux.StepTime.push_back(get_time(t,1,LogAux));
   return ;
}

void read_tracers_ascii(varConfiguration &VarConfigAux)
{
      fprintf(stdout,"\nread_tracers_ascii\n");

   VarConfigAux.NumTrac = 0;
   int NumTot = count_lines(VarConfigAux.FileTracers);
   FILE *fd = safe_open(VarConfigAux.FileTracers,"r");
   float dummy ; // agrego Seba par poder levantar archivo de 7 columnas
   for (int i=0; i<NumTot; i++) {

       //if (NumTrac > 250000) break	   
      fprintf(stdout,"\n Tracer   eliminaer fila\n");

       Tracer.push_back(tracers());

       fscanf(fd,"%f %f %f %f %f %f %f\n",&Tracer[i].Pos[0],&Tracer[i].Pos[1],&Tracer[i].Pos[2],	   
                                        &Tracer[i].Vel[0],&Tracer[i].Vel[1],&Tracer[i].Vel[2],&dummy);
       Tracer[VarConfigAux.NumTrac].Pos[0] *= VarConfigAux.Scale.Pos;
       Tracer[VarConfigAux.NumTrac].Pos[1] *= VarConfigAux.Scale.Pos;
       Tracer[VarConfigAux.NumTrac].Pos[2] *= VarConfigAux.Scale.Pos;
       Tracer[VarConfigAux.NumTrac].Vel[0] *= VarConfigAux.Scale.Vel;
       Tracer[VarConfigAux.NumTrac].Vel[1] *= VarConfigAux.Scale.Vel;
       Tracer[VarConfigAux.NumTrac].Vel[2] *= VarConfigAux.Scale.Vel;
       VarConfigAux.NumTrac++;
   }
   fclose(fd);
   return ;
}

void read_tracers_binary(varConfiguration &VarConfigAux)
{
   FILE *fd = safe_open(VarConfigAux.FileTracers,"r");
   int NumTot;

   VarConfigAux.NumTrac = 0;
   fread(&NumTot,sizeof(int),1,fd);

   for (int i=0; i<NumTot; i++) {
       Tracer.push_back(tracers());

       fread(&Tracer[VarConfigAux.NumTrac].Pos[0],sizeof(float),3,fd);
       fread(&Tracer[VarConfigAux.NumTrac].Vel[0],sizeof(float),3,fd);
       
       Tracer[VarConfigAux.NumTrac].Pos[0] *= VarConfigAux.Scale.Pos;
       Tracer[VarConfigAux.NumTrac].Pos[1] *= VarConfigAux.Scale.Pos;
       Tracer[VarConfigAux.NumTrac].Pos[2] *= VarConfigAux.Scale.Pos;
       Tracer[VarConfigAux.NumTrac].Vel[0] *= VarConfigAux.Scale.Vel;
       Tracer[VarConfigAux.NumTrac].Vel[1] *= VarConfigAux.Scale.Vel;
       Tracer[VarConfigAux.NumTrac].Vel[2] *= VarConfigAux.Scale.Vel;
       VarConfigAux.NumTrac++;
   }
   fclose(fd);
   return ;

}

void read_tracers_gadget(varConfiguration &VarConfigAux)
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
  int    NumFilesAux;

  if (VarConfigAux.NumFiles == 1)
     sprintf(snapshot,"%s",VarConfigAux.FileTracers);
  else 
     sprintf(snapshot,"%s.0",VarConfigAux.FileTracers);

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

  if (Header.BoxSize*VarConfigAux.Scale.Pos != VarConfigAux.BoxSize || Header.NumFiles != VarConfigAux.NumFiles) {
     fprintf(stdout,"\nError. Missmatch with Gadget header.\n");
     fprintf(stdout,"BoxSize = %f (%f in inputfile)\n",Header.BoxSize*VarConfigAux.Scale.Pos,VarConfigAux.BoxSize);
     fprintf(stdout,"NumFiles = %d (%d in inputfile)\n",Header.NumFiles,VarConfigAux.NumFiles);
     fflush(stdout);
     exit(EXIT_FAILURE);
  }

  if ((VarConfigAux.RSDist == 1 || VarConfigAux.GDist == 1) && Header.Redshift != VarConfigAux.Redshift) {
     fprintf(stdout,"\nError. Missmatch with Gadget header.\n");
     fprintf(stdout,"Redshift = %f (%f in inputfile)\n",Header.Redshift,VarConfigAux.Redshift);
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

  VarConfigAux.NumTrac = Header.NpartTotal[1];
  for (int i=0; i<VarConfigAux.NumTrac; i++) Tracer.push_back(tracers());

  if (VarConfigAux.NumFiles < VarConfigAux.OMPcores)
     NC = VarConfigAux.NumFiles;
  else
     NC = VarConfigAux.OMPcores;

  NumFilesAux = VarConfigAux.NumFiles;
  #pragma omp parallel for default(none) schedule(static) num_threads(NC) \
   private(i,snapshot,f1,Np,Header,pos,vel,id,j,k,dummy)  \
   shared(Tracer, VarConfigAux, stdout, NumFilesAux)
  
  for (i=0; i<NumFilesAux; i++) {

      if (NumFilesAux == 1)
         sprintf(snapshot,"%s",VarConfigAux.FileTracers);
      else 
         sprintf(snapshot,"%s.%d",VarConfigAux.FileTracers,i);

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
	      Tracer[id[j]-1].Pos[k] = pos[3*j+k]*VarConfigAux.Scale.Pos;
	      Tracer[id[j]-1].Vel[k] = vel[3*j+k]*sqrt(Header.Time)*VarConfigAux.Scale.Vel;
	  }
      }

      free(pos);
      free(vel);
      free(id);
  }
#undef SKIP 
  return ;
}

void read_tracers_mxxl(varConfiguration &VarConfigAux, logs &LogAux)
{

   int     i,j,k,NumTot,N,id;
   float   Pos[3],Vel[3];
   char    filename[MAXCHAR],basename[MAXCHAR];
   FILE    *fd;
	
   sprintf(basename,"%s/z%4.2f/halos_z%4.2f_part",VarConfigAux.FileTracers,VarConfigAux.Redshift,VarConfigAux.Redshift);
   fprintf(LogAux.logfile," | Files = %s (%d files) \n",basename,VarConfigAux.NumFiles);

   for (i=1; i<=VarConfigAux.NumFiles; i++) {
       sprintf(filename,"%s%02d",basename,i);
       fd = safe_open(filename,"r");
       fread(&N,sizeof(int),1,fd);
       fread(&NumTot,sizeof(int),1,fd);
       VarConfigAux.NumTrac += N;
       fclose(fd);
   }

   if (VarConfigAux.NumTrac != NumTot) {
      fprintf(stdout,"Error. Missing halos: NumTrac = %d, NumTot = %d\n",VarConfigAux.NumTrac,NumTot);
      fflush(stdout);
      exit(EXIT_FAILURE);
   }

   for (int i=0; i<VarConfigAux.NumTrac; i++) Tracer.push_back(tracers());

   for (j=1; j<=VarConfigAux.NumFiles; j++) {
       sprintf(filename,"%s%02d",basename,j);
       fd = safe_open(filename,"r");
       fread(&N,sizeof(int),1,fd);
       fread(&NumTot,sizeof(int),1,fd);
       for (i=0; i<N; i++) {
	   fread(&id,sizeof(int),1,fd);
	   fread(Pos,sizeof(float),3,fd);
	   fread(Vel,sizeof(float),3,fd);

	   for (k=0; k<3; k++) {
	       Tracer[id].Pos[k] = Pos[k]*VarConfigAux.Scale.Pos;
	       Tracer[id].Vel[k] = Vel[k]*VarConfigAux.Scale.Vel;
	   }
       }
       fclose(fd);
   }
   return ;
}

void redshift_space_distortions(varConfiguration VarConfigAux, logs &LogAux)
{

   int    i;
   double Hubble_z;
   double RSDFactor;	
   struct cosmoparam C = {
      VarConfigAux.OmegaMatter,VarConfigAux.OmegaLambda,
      (1.0 - VarConfigAux.OmegaMatter - VarConfigAux.OmegaLambda),
      VarConfigAux.Hubble
   };

   fprintf(LogAux.logfile," | Applying redshift-space distortions: LOS = z-axis, POS = xy-plane\n");

   if (VarConfigAux.Redshift == 0.0) {
     RSDFactor = 1.0/100.0;	   
   } else { 
     Hubble_z = 100.0*evolution_param(VarConfigAux.Redshift,&C);
     RSDFactor = (1.0 + VarConfigAux.Redshift)/Hubble_z;
   }

   fprintf(LogAux.logfile," | RSD factor = (1+z)/H(z) = %f [h⁻¹Mpc/(km/s)]\n",RSDFactor);
   fflush(LogAux.logfile);

   for (i=0; i<VarConfigAux.NumTrac; i++) {
       Tracer[i].Pos[2] += Tracer[i].Vel[2]*RSDFactor;
       if (Tracer[i].Pos[2] < 0.0    ) Tracer[i].Pos[2] += VarConfigAux.LBox[2];
       if (Tracer[i].Pos[2] > VarConfigAux.LBox[2]) Tracer[i].Pos[2] -= VarConfigAux.LBox[2];
   }	       
   return ;
}

void geometrical_distortions(varConfiguration &VarConfigAux, logs &LogAux)
{

   int    i;
   double Hubble_z,Distance_z;	
   double FidHubble_z,FidDistance_z;	
   double GDFactor_LOS,GDFactor_POS;
   struct cosmoparam C = {
      VarConfigAux.OmegaMatter,
      VarConfigAux.OmegaLambda,
		1.0 - VarConfigAux.OmegaMatter - VarConfigAux.OmegaLambda,
      VarConfigAux.Hubble
   };
   struct cosmoparam FC = {
      VarConfigAux.FidOmegaMatter,
	   VarConfigAux.FidOmegaLambda,
		1.0 - VarConfigAux.FidOmegaMatter - VarConfigAux.FidOmegaLambda,
      VarConfigAux.FidHubble
   };

   fprintf(LogAux.logfile," | Applying fiducial cosmology distortions: LOS = z-axis, POS = xy-plane\n");
   fprintf(LogAux.logfile," | True cosmology: (OmM,OmL,OmK,H0) = (%4.2f,%4.2f,%4.2f,%4.2f)\n",C.OmM,C.OmL,C.OmK,C.Hub);

   Hubble_z = C.Hub*evolution_param(VarConfigAux.Redshift,&C);
   Distance_z = angular_distance(VarConfigAux.Redshift,&C);
   
   fprintf(LogAux.logfile," | Fiducial cosmology: (OmM,OmL,OmK,H0) = (%4.2f,%4.2f,%4.2f,%4.2f)\n",FC.OmM,FC.OmL,FC.OmK,FC.Hub);

   FidHubble_z = FC.Hub*evolution_param(VarConfigAux.Redshift,&FC);
   FidDistance_z = angular_distance(VarConfigAux.Redshift,&FC);

   GDFactor_LOS = Hubble_z/FidHubble_z;
   GDFactor_POS = FidDistance_z/Distance_z;

   fprintf(LogAux.logfile," | GD factor: LOS = H(z)/FidH(z) = %f \n",GDFactor_LOS);
   fprintf(LogAux.logfile," |            POS = FidD(z)/D(z) = %f \n",GDFactor_POS);
   fflush(LogAux.logfile);

   VarConfigAux.LBox[0] *= GDFactor_POS;
   VarConfigAux.LBox[1] *= GDFactor_POS;
   VarConfigAux.LBox[2] *= GDFactor_LOS;

   for (i=0; i<VarConfigAux.NumTrac; i++) {
       Tracer[i].Pos[0] *= GDFactor_POS;
       Tracer[i].Pos[1] *= GDFactor_POS;
       Tracer[i].Pos[2] *= GDFactor_LOS;
   }	
   return ;
}

void write_voids(varConfiguration VarConfigAux, logs &LogAux)
{
   int     i;
   FILE    *fd;
   clock_t t;

   fprintf(LogAux.logfile,"\n WRITTING VOID CATALOGUE \n");
   t = clock();
   
   fd = safe_open(VarConfigAux.FileVoids,"w");

   for (i=0; i<VarConfigAux.NumVoid; i++) {
       if (Void[i].ToF) {

          fprintf(fd," %8.5f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %5d\n",
	  	       Void[i].Rad,Void[i].Pos[0],Void[i].Pos[1],Void[i].Pos[2],Void[i].Vel[0],Void[i].Vel[1],
	  	       Void[i].Vel[2],Void[i].Delta,Void[i].Dtype,Void[i].Poisson,Void[i].Nran);   
       }
   }
   fclose(fd);

   LogAux.StepName.push_back("Writting void catalogue");
   LogAux.StepTime.push_back(get_time(t,1,LogAux));

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
