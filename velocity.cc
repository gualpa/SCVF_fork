/**
 * @file allvars.h
 * @brief Calcula las velocidades de expansión de los vacíos cósmicos.
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

#include "allvars.h"
#include "velocity.h"
#include "grid.h"
#include "tools.h"

void compute_velocity(varConfiguration VarConfigAux, logs &LogAux, vector <tracers> TracerAux, vector <voids> &VoidAux)
{
   int          i,k,ic,jc,kc,in,m,NumGrid;
   int          ii,jj,kk,l,next,Counter;
   double       xc[3],xt[3],dx[3],vt[3],GridSize[3];
   double       dist,Radius,PLUS,MaxDist,MinDist,GAP;
   struct grid  *GridList;
   int          NumQuery;
   struct query Query;
   clock_t      t;

   fprintf(LogAux.logfile,"\n COMPUTING VOID BULK VELOCITIES \n");
   t = clock();

   NumGrid = (int)round(cbrt((double)VarConfigAux.NumTrac/100.0));
   if (NumGrid < 50) NumGrid = 50;
   GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid));
   build_grid_list(TracerAux,VarConfigAux.NumTrac,GridList,NumGrid,GridSize,false,VarConfigAux,LogAux);

   // Selecciono grides

   MinDist = 0.0;
   MaxDist = 0.0;
   for (i=0; i<VarConfigAux.NumVoid; i++) {
       if (!VoidAux[i].ToF) continue;
       if (VoidAux[i].Rad > MaxDist) MaxDist = VoidAux[i].Rad;
   }
   MaxDist *= 1.5;
   query_grid(&Query,GridSize,MinDist,MaxDist);
   NumQuery = Query.i.size();
   GAP = 0.5*sqrt(3.0)*max_grid_size(GridSize);
  
   fprintf(LogAux.logfile," | MinDist - MaxDist = %5.3f - %5.3f [Mpc/h], %d grids \n",MinDist,MaxDist,NumQuery);
   fflush(LogAux.logfile);

   #pragma omp parallel for default(none) schedule(dynamic)      \
    shared(VarConfigAux,VoidAux,TracerAux,NumQuery,Query,\
           NumGrid,GridSize,GridList,GAP)          \
   private(i,l,k,m,Radius,xc,ic,jc,kc,ii,jj,kk,next,dx,xt,dist,  \
           Counter,vt,PLUS,in)

   for (i=0; i<VarConfigAux.NumVoid; i++) {
       
       if (!VoidAux[i].ToF) continue;

       Counter = 0;
       PLUS = 0.0;
       Radius = VoidAux[i].Rad;
       for (k=0; k<3; k++) {
           xc[k] = (double)VoidAux[i].Pos[k];
           VoidAux[i].Vel[k] = 0.0;
       }

       ic = (int)(xc[0]/GridSize[0]);
       jc = (int)(xc[1]/GridSize[1]);
       kc = (int)(xc[2]/GridSize[2]);

       do {

          for (in=0; in<NumQuery; in++) {
	    
	      ii = Query.i[in]; 
	      jj = Query.j[in]; 
	      kk = Query.k[in]; 

	      dist = (double)(ii*ii)*(GridSize[0]*GridSize[0])
	           + (double)(jj*jj)*(GridSize[1]*GridSize[1])
	           + (double)(kk*kk)*(GridSize[2]*GridSize[2]);
	      dist = sqrt(dist);

	      if (dist > 1.5*Radius+GAP) continue;

              ii = periodic_grid(ii + ic,NumGrid); 
	      jj = periodic_grid(jj + jc,NumGrid); 
	      kk = periodic_grid(kk + kc,NumGrid); 	  

              l = index_1d(ii,jj,kk,NumGrid);

              if (GridList[l].NumMem == 0) continue;

              for (m=0; m<GridList[l].NumMem; m++) {
		
	          next = GridList[l].Member[m];	

                  for (k=0; k<3; k++) {
                      xt[k] = (double)TracerAux[next].Pos[k];
                      vt[k] = (double)TracerAux[next].Vel[k];
                      dx[k] = periodic_delta(xc[k] - xt[k],VarConfigAux.LBox[k]);
                  }

                  dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
                  dist /= Radius;

                  if (dist > VarConfigAux.InnerShellVel-PLUS && dist < VarConfigAux.OuterShellVel+PLUS) {
                     VoidAux[i].Vel[0] += vt[0];
                     VoidAux[i].Vel[1] += vt[1];
                     VoidAux[i].Vel[2] += vt[2];
                     Counter++;	 
                  }
              } 
          }

          PLUS += 0.05;     

       } while (Counter == 0);

       VoidAux[i].Vel[0] /= (double)Counter;
       VoidAux[i].Vel[1] /= (double)Counter;
       VoidAux[i].Vel[2] /= (double)Counter;
   }
  
   free_grid_list(GridList,NumGrid);
   free_query_grid(&Query);

   LogAux.StepName.push_back("Computing velocities");
   LogAux.StepTime.push_back(get_time(t,VarConfigAux.OMPcores,LogAux));
   return;
}
