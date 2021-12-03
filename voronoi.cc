
#include "allvars.h"
#include "voronoi.h"
#include "tools.h"
#include "grid.h"

void Voronoi()
{

   if (ReadVoronoi == 1) {
   
      fprintf(logfile,"\n READING VORONOI TESSELLATION \n");

      switch (FormatTracers) { 
         
	 case 0:
         ReadVoronoi_ASCII();
         break;
         
      	 case 3:
         ReadVoronoi_MXXL();
         break;

      }

   } else if (ReadVoronoi == 0) {
      
      fprintf(logfile,"\n COMPUTING VORONOI TESSELLATION \n");

      ComputeVoronoi();
      
      if (WriteVoronoi == 1) {
      
	  fprintf(logfile,"\n WRITING VORONOI TESSELLATION \n");

          switch (FormatTracers) {

	     case 0:
	     WriteVoronoi_ASCII();
             break;

             case 3:
             WriteVoronoi_MXXL();
             break;	     

	  }		  

      }      

   }   

}

void ComputeVoronoi()
{
   int            G,i,j,k,p,N,l,id,count,NumGrid;
   struct grid    *GridList;
   vector <int>   IDs,indx;
   double         min[3],max[3],ref[3],GridSize[3];   
   double         xp[3],xc[3],rr,Vol;
   bool           check;
   voronoicell    cell;
   container      *con;
   c_loop_order   *clo;
   particle_order *po;
   FILE           *fd;
   clock_t        t;
  
   t = clock();

   NumGrid = (int)round(cbrt((double)NumTrac/MeanPartPerGrid));
   GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid));
   BuildGridList(Tracer,NumTrac,GridList,NumGrid,GridSize,true);

   Vol = LBox[0]*LBox[1]*LBox[2];

   #pragma omp parallel for default(none) schedule(static)           \
    shared(Vol,stdout,GridSize,NumGrid,LBox,Tracer,GridList,NumTrac) \
    private(l,N,k,j,i,ref,min,max,xp,xc,G,count,p,id,IDs,rr,check,   \
            cell,con,clo,po,indx)                                    \

   for (l=0; l<NumGrid*NumGrid*NumGrid; l++) {

       // Count number of particles
        
       N = 0;
       for (j=0; j<27; j++) {
           i = GridList[l].Neighbour[j];	  
           N += GridList[i].NumMem;
       }
   
       // Global parameters

       G = (int)cbrt((double)N/5.0);

       indx = Index3D(l,NumGrid);
   
       for (k=0; k<3; k++) {
           ref[k] = (double)(indx[k]  )*GridSize[k];
           min[k] = (double)(indx[k]-1)*GridSize[k];
           max[k] = (double)(indx[k]+2)*GridSize[k];	      
       }

       con = new container(min[0],max[0],min[1],max[1],min[2],max[2],G,G,G,false,false,false,8);
       po = new particle_order();
       count = 0;

       // Load particles
       
       for (j=0; j<27; j++) {

           i = GridList[l].Neighbour[j];

           for (p=0; p<GridList[i].NumMem; p++) {
               
               id = GridList[i].Member[p];

     	       for (k=0; k<3; k++) {
     	           xp[k] = (double)Tracer[id].Pos[k];
                   if (xp[k] - ref[k] >  0.5*LBox[k]) xp[k] -= LBox[k];	      
                   if (xp[k] - ref[k] < -0.5*LBox[k]) xp[k] += LBox[k];	      
     	       }

     	       if (i == l) 
                  (*con).put(*po,count,xp[0],xp[1],xp[2]);
     	       else
                  (*con).put(count,xp[0],xp[1],xp[2]);
            
               IDs.push_back(id);
               count++;
           }
       }

       // Create loop class
       
       clo = new c_loop_order(*con,*po);

       // Compute Voronoi cells (only those in *po)

       if ((*clo).start()) do {

          // Id and position of the tracer
          (*clo).pos(i,xp[0],xp[1],xp[2],rr);	 

          check = (*con).compute_cell(cell,*clo);
          if (!check) {
             fprintf(stdout,"\n Error. Voronoi cell can not be computed: i = %d \n",i);
     	     fflush(stdout);
             exit(EXIT_FAILURE);	 
          } 
          
          // Centroid position
          cell.centroid(xc[0],xc[1],xc[2]);

          for (k=0; k<3; k++)  
     	      Tracer[IDs[i]].Cen[k] = (float)PeriodicPos(xp[k] + xc[k],LBox[k]);
         
          // Volume of the cell
          Tracer[IDs[i]].Volume = (float)cell.volume(); 
       
          // Delta
          Tracer[IDs[i]].Delta = (float)((Vol/Tracer[IDs[i]].Volume)/(double)(NumTrac) - 1.0);

       } while ((*clo).inc());

       (*con).clear();
       delete con;
       delete po;
       delete clo;
       IDs.clear();
   }

   FreeGridList(GridList,NumGrid);

   StepName.push_back("Computing Voronoi tessellation");
   StepTime.push_back(Time(t,OMPcores));

}

void ReadVoronoi_ASCII()
{
   FILE    *fd;
   int     i;
   clock_t t;

   t = clock();

   fd = SafeOpen(FileVoronoi,"r");
   for (i=0; i<NumTrac; i++) 
       fscanf(fd,"%f %f %f %f %f \n",&Tracer[i].Cen[0],&Tracer[i].Cen[1],&Tracer[i].Cen[2],
      	                             &Tracer[i].Delta,&Tracer[i].Volume);
   fclose(fd);

   StepName.push_back("Reading Voronoi tessellation");
   StepTime.push_back(Time(t,1));

}

void WriteVoronoi_ASCII()
{

   FILE    *fd;
   int     i;
   clock_t t;

   t = clock();

   fd = SafeOpen(FileVoronoi,"w");
   for (i=0; i<NumTrac; i++)
       fprintf(fd,"%12.6f %12.6f %12.6f %12.6f %12.6f \n",
                  Tracer[i].Cen[0],Tracer[i].Cen[1],Tracer[i].Cen[2],
                  Tracer[i].Delta,Tracer[i].Volume);
   fclose(fd);

   StepName.push_back("Writting Voronoi tessellation");
   StepTime.push_back(Time(t,1));

}

void ReadVoronoi_MXXL()
{
   int     i,j,k,NumFiles,N,NumTot,id;
   float   Vol,Pos[3],Volume;
   FILE    *fd;
   char    filename[MAXCHAR];
   clock_t t;

   Vol = LBox[0]*LBox[1]*LBox[2];

   NumTot = 0;
   for (j=1; j<=NumFiles; j++) { 

       sprintf(filename,"%s.%02d",FileVoronoi,NumFiles);
       fd = SafeOpen(filename,"r");
       fread(&N,sizeof(int),1,fd);

       NumTot += N;

       for (i=0; i<N; i++) {

           fread(&id,sizeof(int),1,fd);
           fread(Pos,sizeof(float),3,fd);
           fread(&Volume,sizeof(float),1,fd);

           for (k=0; k<3; k++) Tracer[id].Cen[k] = (float)Pos[k];
           Tracer[id].Volume = Volume;
           Tracer[id].Delta = (Vol/Tracer[id].Volume)/(float)NumTrac - 1.0;

       }
       fclose(fd);
   }

   if (NumTrac != NumTot) {
      fprintf(stdout,"Error. Missing halos in Voronoi: NumTrac = %d, NumTot = %d\n",NumTrac,NumTot);
      fflush(stdout);
      exit(EXIT_FAILURE);
   }

   StepTime.push_back(Time(t,1));

}

void WriteVoronoi_MXXL()
{

   int     N,i,j,NumTot;
   FILE    *fd;
   char    filename[MAXCHAR];
   clock_t t;

   t = clock();

   N = (int)((double)NumTrac/(double)NumFiles);

   NumTot = 0;

   for (j=1; j<=NumFiles; j++) { 

       sprintf(filename,"%s.%02d",FileVoronoi,NumFiles);
       fd = SafeOpen(filename,"w");

       fwrite(&N,1,sizeof(int),fd);

       for (i=NumTot; i<(NumTot+N); i++) {
	   fwrite(&i,1,sizeof(int),fd);    
	   fwrite(Tracer[i].Cen,3,sizeof(float),fd);    
	   fwrite(&Tracer[i].Volume,1,sizeof(float),fd);    
       }
       fclose(fd);

       NumTot += N;
   }

   if (NumTrac != NumTot) {
      fprintf(stdout,"Error. Missing halos in Voronoi: NumTrac = %d, NumTot = %d\n",NumTrac,NumTot);
      fflush(stdout);
      exit(EXIT_FAILURE);
   }

}

