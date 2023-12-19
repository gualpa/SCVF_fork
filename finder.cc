#include "allvars.h"
#include "finder.h"
#include "grid.h"
#include "qsort.h"
#include "tools.h"

void find_void_candidates()
{
  int     i;
  clock_t t;
   
  fprintf(VarConfig.logfile,"\n SELECTING UNDERDENSE REGIONS \n");
  t = clock();

  VarConfig.NumVoid = 0;
  for (i=0; i<VarConfig.NumTrac; i++) {
      if (Tracer[i].Delta <= VarConfig.DeltaSeed) {
         Void.push_back(voids());
         Void[VarConfig.NumVoid].Pos[0] = Tracer[i].Cen[0];	 
         Void[VarConfig.NumVoid].Pos[1] = Tracer[i].Cen[1];	 
         Void[VarConfig.NumVoid].Pos[2] = Tracer[i].Cen[2];
	 Void[VarConfig.NumVoid].Ini[0] = Tracer[i].Cen[0];	 
         Void[VarConfig.NumVoid].Ini[1] = Tracer[i].Cen[1];	 
         Void[VarConfig.NumVoid].Ini[2] = Tracer[i].Cen[2];
         Void[VarConfig.NumVoid].Rini = 1.5*cbrt(0.75*(double)Tracer[i].Volume/PI);
         Void[VarConfig.NumVoid].Rad = 0.0;
         Void[VarConfig.NumVoid].ToF = false;
	 Void[VarConfig.NumVoid].Nran = 0;
         VarConfig.NumVoid++;	 
      }
  }

  fprintf(VarConfig.logfile," | Void candidates = %d \n",VarConfig.NumVoid);
  
  VarConfig.StepName.push_back("Finding centers");
  VarConfig.StepTime.push_back(get_time(t,1));

}

void find_voids() 
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

  fprintf(VarConfig.logfile,"\n VOID IDENTIFICATION \n");
  t = clock();

  srand(time(NULL));
  
  NumGrid = (int)(VarConfig.BoxSize/VarConfig.ProxyGridSize);
  GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid));
  build_grid_list(Tracer,VarConfig.NumTrac,GridList,NumGrid,GridSize,false);

// int          NumShell = (int)round(0.5*MaxRadiusSearch/max_grid_size(GridSize));
  int          NumShell = (int)round(VarConfig.MaxRadiusSearch/max_grid_size(GridSize));
  int          NumQuery[NumShell];
  struct query Query[NumShell];

  // Selecciono vecinos

  if (VarConfig.OMPcores > NumShell) 
     NumCores = NumShell;
  else
     NumCores = VarConfig.OMPcores;	

  #pragma omp parallel for default(none) num_threads(NumCores)    \
   shared(NumShell,NumQuery,Query,GridSize,stdout,VarConfig.MaxRadiusSearch)\
   private(p,MinDist,MaxDist)

  for (p=0; p<NumShell; p++) {
      MinDist = VarConfig.MaxRadiusSearch/(double)NumShell*(double)p;  
      MaxDist = VarConfig.MaxRadiusSearch/(double)NumShell*(double)(p+1);  
      query_grid(&Query[p],GridSize,MinDist,MaxDist);
      NumQuery[p] = Query[p].i.size();
  } 

  for (p=0; p<NumShell; p++) {
      MinDist = VarConfig.MaxRadiusSearch/(double)NumShell*(double)p;  
      MaxDist = VarConfig.MaxRadiusSearch/(double)NumShell*(double)(p+1); 
      fprintf(VarConfig.logfile," | Shell N° %2d: MinDist - MaxDist = %5.2f - %5.2f [Mpc/h], %5d grids (Overlap = %f) \n",
		      p,MinDist,MaxDist,NumQuery[p],0.5*sqrt(3.0)*max_grid_size(GridSize));
  }
  fflush(VarConfig.logfile);

  #pragma omp parallel for default(none) schedule(static)                     \
   shared(VarConfig.NumVoid,VarConfig.MeanNumTrac,Void,Tracer,Query,NumQuery,NumShell,VarConfig.LBox,stdout,\
          VarConfig.MaxRadiusSearch,VarConfig.FracRadius,VarConfig.DeltaThreshold,NumGrid,GridSize,GridList,\
	         VarConfig.OMPcores,VarConfig.RadIncrement,NumRanWalk)                                   \
   private(iv,ir,ic,jc,kc,ii,jj,kk,xc,xr,dx,l,Radius,BiggestRadius,next,k,    \
           dist,val,kappa,SortArr,Nsort,the,phi,rad,Volume,Delta,lambda,p,m,  \
	   MinDist,MaxDist,done,in,CheckRan,TotRan)

  for (iv=0; iv<VarConfig.NumVoid; iv++) {

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
	     rad = (double)Void[iv].Rad;

	     if (rad == 0.0) 
	        rad = (double)Void[iv].Rini*random_number();
	     else  
	        rad *= VarConfig.FracRadius*random_number();

	     xr[0] = rad*sin(the)*cos(phi);
	     xr[1] = rad*sin(the)*sin(phi);
	     xr[2] = rad*cos(the);
	  }

	  if (Void[iv].Rad == 0.0) {

             for (k=0; k<3; k++) 
		 xc[k] = periodic_position((double)Void[iv].Ini[k] + xr[k],VarConfig.LBox[k]);

	  } else {

             for (k=0; k<3; k++) {
	         xc[k] = periodic_position((double)Void[iv].Pos[k] + xr[k],VarConfig.LBox[k]);
	         dx[k] = periodic_delta(xc[k] - (double)Void[iv].Ini[k],VarConfig.LBox[k]);
	     }
	     dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
	     if (dist > Void[iv].Rad) // avoid big migration 
	        for (k=0; k<3; k++) 
	            xc[k] = periodic_position((double)Void[iv].Ini[k] + xr[k],VarConfig.LBox[k]); 	  
	  }

	  ic = (int)(xc[0]/GridSize[0]);
	  jc = (int)(xc[1]/GridSize[1]);
	  kc = (int)(xc[2]/GridSize[2]);

	  done = false;
	  p = 0;

	  do {

             MinDist = VarConfig.MaxRadiusSearch/(double)NumShell*(double)p;		  
             MaxDist = VarConfig.MaxRadiusSearch/(double)NumShell*(double)(p+1);		  

	     for (in=0; in<NumQuery[p]; in++) {

	         ii = periodic_grid(Query[p].i[in] + ic,NumGrid); 
	         jj = periodic_grid(Query[p].j[in] + jc,NumGrid); 
	         kk = periodic_grid(Query[p].k[in] + kc,NumGrid); 

		 l = index_1d(ii,jj,kk,NumGrid);

		 if (GridList[l].NumMem == 0) continue;

		 for (m=0; m<GridList[l].NumMem; m++) {
		
	             next = GridList[l].Member[m];		 

		     for (k=0; k<3; k++) 
	                 dx[k] = periodic_delta(xc[k] - (double)Tracer[next].Pos[k],VarConfig.LBox[k]);

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
	     Delta = (double)(Nsort-1)/Volume/VarConfig.MeanNumTrac - 1.0;

	     if (Delta < VarConfig.DeltaThreshold) {
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
	      Delta = (double)(ir+1)/Volume/VarConfig.MeanNumTrac - 1.0;

              if (Delta < VarConfig.DeltaThreshold && Radius > BiggestRadius) {
		 
		 if (Radius/BiggestRadius - 1.0 >= VarConfig.RadIncrement) CheckRan = 0;
		
	         Void[iv].Rad = (float)Radius;
		 Void[iv].Delta = (float)Delta;
	         Void[iv].Pos[0] = (float)xc[0];	       
	         Void[iv].Pos[1] = (float)xc[1];	       
	         Void[iv].Pos[2] = (float)xc[2];
	         Void[iv].ToF = true;
		 kappa = ir + 1; 
	         BiggestRadius = Radius;
	      
	      } /* Fin lazo Dcum < VarConfig.DeltaThreshold */

	  } /* Fin lazo bines */
	 
	  free(SortArr);
 
      } while (CheckRan < VarConfig.NumRanWalk); /* Fin lazo random */      

      if (Void[iv].ToF) {
	 lambda = (4.0/3.0)*PI*pow((double)Void[iv].Rad,3)*VarConfig.MeanNumTrac;
         Void[iv].Poisson = (double)kappa*log(lambda) - lambda - ln_factorial(kappa); 	 
	 Void[iv].Nran = TotRan;
      }    

  } /* Fin lazo voids */

  for (p=0; p<NumShell; p++) 
      free_query_grid(&Query[p]);	  
  free_grid_list(GridList,NumGrid);

  VarConfig.StepName.push_back("Finding voids");
  VarConfig.StepTime.push_back(get_time(t,VarConfig.OMPcores));

}

void clean_voids()
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

  fprintf(VarConfig.logfile,"\n CLEANING VOID CATALOGUE BY OVERLAP (TOL = %4.2f) \n",OverlapTol);
  t = clock();

  NumGrid = (int)cbrt((double)VarConfig.NumVoid/10.0);
  GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid)); 
  build_grid_list(Void,VarConfig.NumVoid,GridList,NumGrid,GridSize,false);
  
  NumTrueVoid = 0;
  for (i=0; i<VarConfig.NumVoid; i++)
      if (Void[i].ToF) 
	 NumTrueVoid++;	  
  
  fprintf(VarConfig.logfile," | Number of true voids before cleaning = %d (%4.2f %)\n",
		  NumTrueVoid,(double)NumTrueVoid/(double)VarConfig.NumVoid*100.0);
  fflush(VarConfig.logfile);

  SortArr = (struct sort *) malloc(NumTrueVoid*sizeof(struct sort));

  it = 0;
  for (i=0; i<VarConfig.NumVoid; i++) {
      if (Void[i].ToF) {	      
         SortArr[it].val = Void[i].Rad;
         SortArr[it].ord = i;
	 it++;
      }     
  }

  qsort(SortArr,0,NumTrueVoid-1);

  // Selecciono vecinos
  
  MinDist = 0.0;
  MaxDist = 2.0*VarConfig.MaxRadiusSearch;  

  query_grid(&Query,GridSize,MinDist,MaxDist);
  NumQuery = Query.i.size();
  
  fprintf(VarConfig.logfile," | MinDist - MaxDist = %5.3f - %5.3f [Mpc/h], %d grids \n",MinDist,MaxDist,NumQuery);
  fflush(VarConfig.logfile);

  for (i=NumTrueVoid-1; i>=0; i--) {

      indx = SortArr[i].ord;
      if (!Void[indx].ToF) continue;

      xi[0] = (double)Void[indx].Pos[0];      
      xi[1] = (double)Void[indx].Pos[1];      
      xi[2] = (double)Void[indx].Pos[2];
      Ri = (double)Void[indx].Rad;
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

	      if (Void[next].ToF && next != indx) {

                 xj[0] = (double)Void[next].Pos[0];      
                 xj[1] = (double)Void[next].Pos[1];      
                 xj[2] = (double)Void[next].Pos[2];
                 Rj = (double)Void[next].Rad;
	         Vj = (4.0/3.0)*PI*Rj*Rj*Rj;
    
	         for (k=0; k<3; k++) 
	             dx[k] = periodic_delta(xi[k] - xj[k],VarConfig.LBox[k]);

                 dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

	         // Volumen de la interseccion, normalizado a Vj
                 // Solo vale para d <= Ri+Rj

	         if (dist > Ri + Rj) continue;

	         Vij  = (PI/12.0/dist);
                 Vij *= pow(Ri + Rj - dist,2);
	         Vij *= (pow(dist,2) + 2.0*dist*(Ri + Rj) - 3.0*pow(Ri - Rj,2));
	         Vij /= Vj;

	         if (Vij > VarConfig.OverlapTol) Void[next].ToF = false;
    
	      } 
    	  }
      }
  }

  NumTrueVoid = 0;
  for (i=0; i<VarConfig.NumVoid; i++)
      if (Void[i].ToF) 
	 NumTrueVoid++;	  
  
  fprintf(VarConfig.logfile," | Number of true void after cleaning = %d (%4.2f %)\n",
		  NumTrueVoid,(double)NumTrueVoid/(double)VarConfig.NumVoid*100.0);

  free(SortArr);
  free_grid_list(GridList,NumGrid);
  free_query_grid(&Query);

  VarConfig.StepName.push_back("Cleaning voids"); 
  VarConfig.StepTime.push_back(get_time(t,1));

}

