
#include "allvars.h"
#include "proto.h"

void ComputeVelocity()
{

   int              i,j,k,ic,jc,kc,in,m,NumGrid;
   int              ii,jj,kk,l,next,Counter,NG[3];
   double           xc[3],xt[3],dx[3],vt[3],GAP,GridSize[3];
   double           dist,Radius,PLUS,MaxDist,MinDist;
   struct grid      *GridList;
   int              NumNeigh;
   struct neighbour Neigh;
   clock_t          t;

   fprintf(logfile,"\n COMPUTING VOID BULK VELOCITIES \n");
   t = clock();

   NumGrid = (int)round(cbrt((double)NumTrac/100.0));
   if (NumGrid < 50) NumGrid = 50;
   GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid));
   BuildGridList(GridList,NumGrid,GridSize,0,false);

   GAP = 0.0;
   for (k=0; k<3; k++) 
       if (GridSize[k] > GAP) 
          GAP = GridSize[k];	  
   GAP *= sqrt(3.0);

   MinDist = 0.0;
   MaxDist = OuterShellVel*MaxRadiusSearch + GAP;  	 

   SearchNeighbours(&Neigh,&NumNeigh,GridSize,MinDist,MaxDist);
  
   fprintf(logfile," | MinDist - MaxDist = %5.3f - %5.3f [Mpc/h], %d grids \n",MinDist,MaxDist,NumNeigh);
   fflush(logfile);

   #pragma omp parallel for default(none) schedule(dynamic)                    \
    shared(NumVoid,Void,Tracer,NumNeigh,Neigh,LBox,InnerShellVel,OuterShellVel,\
           NumGrid,GridSize,GridList)                                          \
   private(i,l,k,m,Radius,xc,ic,jc,kc,ii,jj,kk,next,dx,xt,dist,Counter,vt,PLUS,in)

   for (i=0; i<NumVoid; i++) {
       
       //if (omp_get_thread_num() == 0) Progress(i,NumVoid);
       
       if (!Void[i].ToF) continue;

       Counter = 0;
       PLUS = 0.0;
       Radius = Void[i].Rad;
       for (k=0; k<3; k++) {
           xc[k] = (double)Void[i].Pos[k];
           Void[i].Vel[k] = 0.0;
       }

       ic = (int)(xc[0]/GridSize[0]);
       jc = (int)(xc[1]/GridSize[1]);
       kc = (int)(xc[2]/GridSize[2]);

       do {

          for (in=0; in<NumNeigh; in++) {
	    
              ii = PeriodicGrid(Neigh.i[in] + ic,NumGrid); 
	      jj = PeriodicGrid(Neigh.j[in] + jc,NumGrid); 
	      kk = PeriodicGrid(Neigh.k[in] + kc,NumGrid); 	  

              l = Index1D(ii,jj,kk,NumGrid);

              if (GridList[l].NumMem == 0) continue;

              for (m=0; m<GridList[l].NumMem; m++) {
		
	          next = GridList[l].Member[m];	

                  for (k=0; k<3; k++) {
                      xt[k] = (double)Tracer[next].Pos[k];	 
                      vt[k] = (double)Tracer[next].Vel[k];
                      dx[k] = PeriodicDeltaPos(xc[k] - xt[k],LBox[k]);
                  }

                  dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
                  dist /= Radius;

                  if (dist > InnerShellVel-PLUS && dist < OuterShellVel+PLUS) {
                     Void[i].Vel[0] += vt[0];
                     Void[i].Vel[1] += vt[1];
                     Void[i].Vel[2] += vt[2];
                     Counter++;	 
                  }
              } 
          }

          PLUS += 0.05;     

       } while (Counter == 0);

       Void[i].Vel[0] /= (double)Counter;
       Void[i].Vel[1] /= (double)Counter;
       Void[i].Vel[2] /= (double)Counter; 
   }
  
   FreeGridList(GridList,NumGrid);
   FreeNeighbours(&Neigh);

   StepName.push_back("Computing velocities");
   StepTime.push_back(Time(t,OMPcores));
}
