/**
 * @file allvars.h
 * @brief Este archivo contiene funciones para identificar y caracterizar los "voids" en
 *        una distribución espacial 3D de trazadores.
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
#include "finder.h"
#include "grid.h"
#include "qsort.h"
#include "tools.h"

vector <voids> find_void_candidates(varConfiguration &VarConfigAux, logs &LogAux, vector <tracers> TracerAux)
{
  int     i;
  clock_t t;
  vector <voids> VoidAux;

  fprintf(LogAux.logfile,"\n SELECTING UNDERDENSE REGIONS \n");
  t = clock();

  VarConfigAux.NumVoid = 0;
  for (i=0; i<VarConfigAux.NumTrac; i++) {

      if (TracerAux[i].Delta <= VarConfigAux.DeltaSeed) {

         VoidAux.push_back(voids());
         VoidAux[VarConfigAux.NumVoid].Pos[0] = TracerAux[i].Cen[0];
         VoidAux[VarConfigAux.NumVoid].Pos[1] = TracerAux[i].Cen[1];
         VoidAux[VarConfigAux.NumVoid].Pos[2] = TracerAux[i].Cen[2];
	       VoidAux[VarConfigAux.NumVoid].Ini[0] = TracerAux[i].Cen[0];
         VoidAux[VarConfigAux.NumVoid].Ini[1] = TracerAux[i].Cen[1];
         VoidAux[VarConfigAux.NumVoid].Ini[2] = TracerAux[i].Cen[2];
         VoidAux[VarConfigAux.NumVoid].Rini = 1.5*cbrt(0.75*(double)TracerAux[i].Volume/PI);
         VoidAux[VarConfigAux.NumVoid].Rad = 0.0;
         VoidAux[VarConfigAux.NumVoid].ToF = false;
	       VoidAux[VarConfigAux.NumVoid].Nran = 0;
         VarConfigAux.NumVoid++;
      }
  }

  fprintf(LogAux.logfile," | Void candidates = %d \n",VarConfigAux.NumVoid);
  
  LogAux.StepName.push_back("Finding centers");
  LogAux.StepTime.push_back(get_time(t,1,LogAux));
  return VoidAux;
}

vector <voids> find_voids(varConfiguration VarConfigAux, logs &LogAux, vector <tracers> TracerAux, vector <voids> VoidAux)
{
  struct grid    *GridList;
  int            NumCores,NumGrid;
  int            iv,CheckRan,TotRan,ir,next,l,kappa,p,in;
  int            ic,jc,kc,ii,jj,kk,k,Nsort,m;
  double         Radius,BiggestRadius,lambda,MinDist,MaxDist;
  double         dx[3],xr[3],xc[3],dist,GridSize[3];
  double         the,phi,rad,Volume,Delta;
  vector <float> val;
  struct sort    *SortArr;
  bool           done;
  clock_t        t;

  fprintf(LogAux.logfile,"\n VOID IDENTIFICATION \n");
  t = clock();
  srand(time(NULL));

  NumGrid = (int)(VarConfigAux.BoxSize/VarConfigAux.ProxyGridSize);
  GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid));
  build_grid_list(TracerAux,VarConfigAux.NumTrac,GridList,NumGrid,GridSize,false,VarConfigAux,LogAux);
// int          NumShell = (int)round(0.5*MaxRadiusSearch/max_grid_size(GridSize));
  int          NumShell = (int)round(VarConfigAux.MaxRadiusSearch/max_grid_size(GridSize));
  int          NumQuery[NumShell];
  struct query Query[NumShell];
  // Selecciono vecinos
  if (VarConfigAux.OMPcores > NumShell)
     NumCores = NumShell;
  else
     NumCores = VarConfigAux.OMPcores;
  #pragma omp parallel for default(none) num_threads(NumCores)    \
   shared(NumShell,NumQuery,Query,GridSize,stdout,VarConfigAux,LogAux)\
   private(p,MinDist,MaxDist)
  for (p=0; p<NumShell; p++) {
      MinDist = VarConfigAux.MaxRadiusSearch/(double)NumShell*(double)p;
      MaxDist = VarConfigAux.MaxRadiusSearch/(double)NumShell*(double)(p+1);
      query_grid(&Query[p],GridSize,MinDist,MaxDist);
      NumQuery[p] = Query[p].i.size();
      fprintf(LogAux.logfile," | Shell N° %2d: MinDist - MaxDist = %5.2f - %5.2f [Mpc/h], %5d grids (Overlap = %f) \n",
              p,MinDist,MaxDist,NumQuery[p],0.5*sqrt(3.0)*max_grid_size(GridSize));
  } 

  /*for (p=0; p<NumShell; p++) {
      MinDist = VarConfigAux.MaxRadiusSearch/(double)NumShell*(double)p;
      MaxDist = VarConfigAux.MaxRadiusSearch/(double)NumShell*(double)(p+1);
      fprintf(LogAux.logfile," | Shell N° %2d: MinDist - MaxDist = %5.2f - %5.2f [Mpc/h], %5d grids (Overlap = %f) \n",
		      p,MinDist,MaxDist,NumQuery[p],0.5*sqrt(3.0)*max_grid_size(GridSize));
  }*/
  fflush(LogAux.logfile);


  #pragma omp parallel for default(none) schedule(static)                     \
   shared(VarConfigAux,VoidAux,TracerAux,Query,NumQuery,NumShell,stdout,\
          NumGrid,GridSize,GridList)                                   \
   private(iv,ir,ic,jc,kc,ii,jj,kk,xc,xr,dx,l,Radius,BiggestRadius,next,k,    \
           dist,val,kappa,SortArr,Nsort,the,phi,rad,Volume,Delta,lambda,p,m,  \
	   MinDist,MaxDist,done,in,CheckRan,TotRan)

  for (iv=0; iv<VarConfigAux.NumVoid; iv++) {

      BiggestRadius = 0.1; 
      TotRan = 0;
      CheckRan = 0;

      do {

	  TotRan++;
	  CheckRan++;

	  if (TotRan == 1) {

             for (k=0; k<3; k++) 
		 xr[k] = 0.0; 

	  } else {

	     the = acos(2.0*random_number() - 1.0);
	     phi = 2.0*PI*random_number();
	     rad = (double)VoidAux[iv].Rad;

	     if (rad == 0.0) 
	        rad = (double)VoidAux[iv].Rini*random_number();
	     else  
	        rad *= VarConfigAux.FracRadius*random_number();

	     xr[0] = rad*sin(the)*cos(phi);
	     xr[1] = rad*sin(the)*sin(phi);
	     xr[2] = rad*cos(the);
	  }

	  if (VoidAux[iv].Rad == 0.0) {

             for (k=0; k<3; k++) 
		 xc[k] = periodic_position((double)VoidAux[iv].Ini[k] + xr[k],VarConfigAux.LBox[k]);

	  } else {

             for (k=0; k<3; k++) {
	         xc[k] = periodic_position((double)VoidAux[iv].Pos[k] + xr[k],VarConfigAux.LBox[k]);
	         dx[k] = periodic_delta(xc[k] - (double)VoidAux[iv].Ini[k],VarConfigAux.LBox[k]);
	     }
	     dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
	     if (dist > VoidAux[iv].Rad) // avoid big migration
	        for (k=0; k<3; k++) 
	            xc[k] = periodic_position((double)VoidAux[iv].Ini[k] + xr[k],VarConfigAux.LBox[k]);
	  }

	  ic = (int)(xc[0]/GridSize[0]);
	  jc = (int)(xc[1]/GridSize[1]);
	  kc = (int)(xc[2]/GridSize[2]);

	  done = false;
	  p = 0;

	  do {

             MinDist = VarConfigAux.MaxRadiusSearch/(double)NumShell*(double)p;
             MaxDist = VarConfigAux.MaxRadiusSearch/(double)NumShell*(double)(p+1);

	     for (in=0; in<NumQuery[p]; in++) {

	         ii = periodic_grid(Query[p].i[in] + ic,NumGrid); 
	         jj = periodic_grid(Query[p].j[in] + jc,NumGrid); 
	         kk = periodic_grid(Query[p].k[in] + kc,NumGrid); 

		 l = index_1d(ii,jj,kk,NumGrid);

		 if (GridList[l].NumMem == 0) continue;

		 for (m=0; m<GridList[l].NumMem; m++) {
		
	             next = GridList[l].Member[m];		 

		     for (k=0; k<3; k++) 
	                 dx[k] = periodic_delta(xc[k] - (double)TracerAux[next].Pos[k],VarConfigAux.LBox[k]);

                     dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

	             if (dist >= MinDist && dist < MaxDist)
	                val.push_back(dist);

		 }

	     }

	     Nsort = val.size();
	     SortArr = (struct sort *) malloc(Nsort*sizeof(struct sort));
	     for (k=0; k<Nsort; k++) {
	         SortArr[k].val = val[k];
                 SortArr[k].ord = k;	      
	     }
             qsort(SortArr,0,Nsort-1);

	     Radius = 0.5*(SortArr[Nsort-2].val + SortArr[Nsort-1].val);    
   	     Volume = (4.0/3.0)*PI*Radius*Radius*Radius;
	     Delta = (double)(Nsort-1)/Volume/VarConfigAux.MeanNumTrac - 1.0;

	     if (Delta < VarConfigAux.DeltaThreshold) {
	        p++;    
		free(SortArr); 
	     } else {
		done = true;
		val.clear();
	     }            

	  } while(!done && p != NumShell);

	  if (!done && p == NumShell) {
             fprintf(stdout,"\n Error. Increase MaxRadiusSearch.\n");
	     fflush(stdout);
             exit(EXIT_FAILURE);	     
	  }

	  for (ir=0; ir<Nsort-1; ir++) {
	
	      Radius = 0.5*(SortArr[ir].val + SortArr[ir+1].val);    
   	      Volume = (4.0/3.0)*PI*Radius*Radius*Radius;
	      Delta = (double)(ir+1)/Volume/VarConfigAux.MeanNumTrac - 1.0;

        if (Delta < VarConfigAux.DeltaThreshold && Radius > BiggestRadius) {
		 
		      if (Radius/BiggestRadius - 1.0 >= VarConfigAux.RadIncrement) CheckRan = 0;
		
          VoidAux[iv].Rad = (float)Radius;
          VoidAux[iv].Delta = (float)Delta;
	        VoidAux[iv].Pos[0] = (float)xc[0];
	        VoidAux[iv].Pos[1] = (float)xc[1];
	        VoidAux[iv].Pos[2] = (float)xc[2];
	        VoidAux[iv].ToF = true;
		      kappa = ir + 1;
	        BiggestRadius = Radius;
	      
	      } /* Fin lazo Dcum < VarConfigAux.DeltaThreshold */

	    } /* Fin lazo bines */
	 
	      free(SortArr);
 
      } while (CheckRan < VarConfigAux.NumRanWalk); /* Fin lazo random */

      if (VoidAux[iv].ToF) {
	 lambda = (4.0/3.0)*PI*pow((double)VoidAux[iv].Rad,3)*VarConfigAux.MeanNumTrac;
         VoidAux[iv].Poisson = (double)kappa*log(lambda) - lambda - ln_factorial(kappa);
	 VoidAux[iv].Nran = TotRan;
      }    

  } /* Fin lazo voids */

  for (p=0; p<NumShell; p++) 
      free_query_grid(&Query[p]);	  
  free_grid_list(GridList,NumGrid);

  LogAux.StepName.push_back("Finding voids");
  LogAux.StepTime.push_back(get_time(t,VarConfigAux.OMPcores,LogAux));

  return VoidAux;
}

vector <voids> clean_voids(varConfiguration VarConfigAux, logs &LogAux, vector <voids> &VoidAux)
{
  int          i,k,l,m,indx,next,in,it;
  int          ii,jj,kk,ic,jc,kc,NumTrueVoid;
  int          NumQuery,NumGrid;
  double       dist,Ri,Rj,Vi,Vj,Vij,MinDist,MaxDist;
  double       xi[3],xj[3],dx[3],GridSize[3];
  struct query Query;
  struct sort  *SortArr; 
  struct grid  *GridList;
  clock_t      t;

  fprintf(LogAux.logfile,"\n CLEANING VOID CATALOGUE BY OVERLAP (TOL = %4.2f) \n",VarConfigAux.OverlapTol);
  t = clock();

  NumGrid = (int)cbrt((double)VarConfigAux.NumVoid/10.0);
  GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid)); 
  build_grid_list(VoidAux,VarConfigAux.NumVoid,GridList,NumGrid,GridSize,false,VarConfigAux,LogAux);
  
  NumTrueVoid = 0;
  for (i=0; i<VarConfigAux.NumVoid; i++)
      if (VoidAux[i].ToF)
	 NumTrueVoid++;	  
  
  fprintf(LogAux.logfile," | Number of true voids before cleaning = %d (%4.2f %)\n",
		  NumTrueVoid,(double)NumTrueVoid/(double)VarConfigAux.NumVoid*100.0);
  fflush(LogAux.logfile);

  SortArr = (struct sort *) malloc(NumTrueVoid*sizeof(struct sort));

  it = 0;
  for (i=0; i<VarConfigAux.NumVoid; i++) {
      if (VoidAux[i].ToF) {
         SortArr[it].val = VoidAux[i].Rad;
         SortArr[it].ord = i;
	 it++;
      }     
  }

  qsort(SortArr,0,NumTrueVoid-1);

  // Selecciono vecinos
  
  MinDist = 0.0;
  MaxDist = 2.0*VarConfigAux.MaxRadiusSearch;

  query_grid(&Query,GridSize,MinDist,MaxDist);
  NumQuery = Query.i.size();
  
  fprintf(LogAux.logfile," | MinDist - MaxDist = %5.3f - %5.3f [Mpc/h], %d grids \n",MinDist,MaxDist,NumQuery);
  fflush(LogAux.logfile);

  for (i=NumTrueVoid-1; i>=0; i--) {

      indx = SortArr[i].ord;
      if (!VoidAux[indx].ToF) continue;

      xi[0] = (double)VoidAux[indx].Pos[0];
      xi[1] = (double)VoidAux[indx].Pos[1];
      xi[2] = (double)VoidAux[indx].Pos[2];
      Ri = (double)VoidAux[indx].Rad;
      Vi = (4.0/3.0)*PI*Ri*Ri*Ri;

      ic = (int)(xi[0]/GridSize[0]);
      jc = (int)(xi[1]/GridSize[1]);
      kc = (int)(xi[2]/GridSize[2]);
    
      for (in=0; in<NumQuery; in++) {

   	  ii = periodic_grid(Query.i[in] + ic,NumGrid); 
	  jj = periodic_grid(Query.j[in] + jc,NumGrid); 
	  kk = periodic_grid(Query.k[in] + kc,NumGrid); 

	  l = index_1d(ii,jj,kk,NumGrid);

          if (GridList[l].NumMem == 0) continue;

	  for (m=0; m<GridList[l].NumMem; m++) {
		
	      next = GridList[l].Member[m];		

	      if (VoidAux[next].ToF && next != indx) {

                 xj[0] = (double)VoidAux[next].Pos[0];
                 xj[1] = (double)VoidAux[next].Pos[1];
                 xj[2] = (double)VoidAux[next].Pos[2];
                 Rj = (double)VoidAux[next].Rad;
	         Vj = (4.0/3.0)*PI*Rj*Rj*Rj;
    
	         for (k=0; k<3; k++) 
	             dx[k] = periodic_delta(xi[k] - xj[k],VarConfigAux.LBox[k]);

                 dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

	         // Volumen de la interseccion, normalizado a Vj
                 // Solo vale para d <= Ri+Rj

	         if (dist > Ri + Rj) continue;

	         Vij  = (PI/12.0/dist);
                 Vij *= pow(Ri + Rj - dist,2);
	         Vij *= (pow(dist,2) + 2.0*dist*(Ri + Rj) - 3.0*pow(Ri - Rj,2));
	         Vij /= Vj;

	         if (Vij > VarConfigAux.OverlapTol) VoidAux[next].ToF = false;
    
	      } 
    	  }
      }
  }

  NumTrueVoid = 0;
  for (i=0; i<VarConfigAux.NumVoid; i++)
      if (VoidAux[i].ToF)
	 NumTrueVoid++;	  
  
  fprintf(LogAux.logfile," | Number of true void after cleaning = %d (%4.2f %)\n",
		  NumTrueVoid,(double)NumTrueVoid/(double)VarConfigAux.NumVoid*100.0);

  free(SortArr);
  free_grid_list(GridList,NumGrid);
  free_query_grid(&Query);

  LogAux.StepName.push_back("Cleaning voids");
  LogAux.StepTime.push_back(get_time(t,1,LogAux));
  return VoidAux;
}

