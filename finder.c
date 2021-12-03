#include "allvars.h"
#include "proto.h"

void FindCenters()
{
  int     i;
  clock_t t;
   
  fprintf(logfile,"\n SELECTING UNDERDENSE REGIONS \n");
  t = clock();

  NumVoid = 0;
  for (i=0; i<NumTrac; i++) {
      if (Tracer[i].Delta <= DeltaSeed) {
         Void.push_back(voids());
         Void[NumVoid].Pos[0] = Tracer[i].Cen[0];	 
         Void[NumVoid].Pos[1] = Tracer[i].Cen[1];	 
         Void[NumVoid].Pos[2] = Tracer[i].Cen[2];
	 Void[NumVoid].Ini[0] = Tracer[i].Cen[0];	 
         Void[NumVoid].Ini[1] = Tracer[i].Cen[1];	 
         Void[NumVoid].Ini[2] = Tracer[i].Cen[2];
         Void[NumVoid].Rini = 1.5*cbrt(0.75*(double)Tracer[i].Volume/PI);
         Void[NumVoid].Rad = 0.0;
         Void[NumVoid].ToF = false;
	 Void[NumVoid].Nran = 0;
         NumVoid++;	 
      }
  }

  fprintf(logfile," | Void candidates = %d \n",NumVoid);
  
  StepName.push_back("Finding centers");
  StepTime.push_back(Time(t,1));

}

void FindVoids() 
{
  struct grid    *GridList;
  int            NumCores,NumGrid;
  int            iv,CheckRan,TotRan,ir,next,l,NC,kappa,p,in;
  int            ic,jc,kc,ii,jj,kk,i,j,k,Nsort,m,NN;
  double         Radius,BiggestRadius,lambda,MinDist,MaxDist;
  double         dx[3],xr[3],xc[3],dist,x,y,z,GridSize[3];
  double         the,phi,rad,Volume,Delta,GAP;
  vector <float> val;
  struct sort    *SortArr;
  bool           done;
  clock_t        t;

  fprintf(logfile,"\n VOID IDENTIFICATION \n");
  t = clock();

  srand(time(NULL));
  
  NumGrid = (int)round(cbrt((double)NumTrac/10.0));
  if (NumGrid < 100) NumGrid = 100;
  GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid));
  BuildGridList(GridList,NumGrid,GridSize,0,false);

  GAP = 0.0;
  for (k=0; k<3; k++) 
      if (GridSize[k] > GAP) 
	 GAP = GridSize[k];	  
  GAP *= sqrt(3.0);
  
  int              NumShell = (int)round(MaxRadiusSearch/0.5/GAP);
  int              NumNeigh[NumShell];
  struct neighbour Neigh[NumShell];

  // Selecciono vecinos

  if (OMPcores > NumShell) 
     NumCores = NumShell;
  else
     NumCores = OMPcores;	

  #pragma omp parallel for default(none) num_threads(NumCores)        \
   shared(NumShell,NumNeigh,Neigh,GridSize,stdout,GAP,MaxRadiusSearch)\
   private(p,MinDist,MaxDist,NN)

  for (p=0; p<NumShell; p++) {

      MinDist = MaxRadiusSearch/(double)NumShell*(double)(p  ) - GAP;  
      MaxDist = MaxRadiusSearch/(double)NumShell*(double)(p+1) + GAP;  

      SearchNeighbours(&Neigh[p],&NumNeigh[p],GridSize,MinDist,MaxDist);
  } 

  for (p=0; p<NumShell; p++) {
      MinDist = MaxRadiusSearch/(double)NumShell*(double)(p  ) - GAP;  
      MaxDist = MaxRadiusSearch/(double)NumShell*(double)(p+1) + GAP; 
      fprintf(logfile," | Shell N° %2d: MinDist - MaxDist = %5.2f - %5.2f [Mpc/h], %5d grids \n",
		      p,MinDist,MaxDist,NumNeigh[p]);
  }
  fflush(logfile);

  #pragma omp parallel for default(none) schedule(static)                     \
   shared(NumVoid,MeanNumTrac,Void,Tracer,Neigh,NumNeigh,NumShell,LBox,stdout,\
          MaxRadiusSearch,FracRadius,DeltaThreshold,NumGrid,GridSize,GridList,\
	  OMPcores,RadIncrement,NumRanWalk,GAP)                               \
   private(iv,ir,ic,jc,kc,ii,jj,kk,xc,xr,dx,l,Radius,BiggestRadius,next,k,    \
           dist,val,kappa,SortArr,Nsort,the,phi,rad,Volume,Delta,lambda,p,m,  \
	   MinDist,MaxDist,done,in,CheckRan,TotRan)

  for (iv=0; iv<NumVoid; iv++) {

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

	     the = acos(2.0*RandomNumber() - 1.0);
	     phi = 2.0*PI*RandomNumber();
	     rad = (double)Void[iv].Rad;

	     if (rad == 0.0) 
	        rad = (double)Void[iv].Rini*RandomNumber();
	     else  
	        rad *= FracRadius*RandomNumber();

	     xr[0] = rad*sin(the)*cos(phi);
	     xr[1] = rad*sin(the)*sin(phi);
	     xr[2] = rad*cos(the);
	  }

	  if (Void[iv].Rad == 0.0) {

             for (k=0; k<3; k++) 
		 xc[k] = PeriodicPos((double)Void[iv].Ini[k] + xr[k],LBox[k]);

	  } else {

             for (k=0; k<3; k++) {
		 xc[k] = PeriodicPos((double)Void[iv].Pos[k] + xr[k],LBox[k]);
		 dx[k] = PeriodicDeltaPos(xc[k] - (double)Void[iv].Ini[k],LBox[k]);
	     }
	     dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
	     if (dist > Void[iv].Rad) // avoid big migration 
	        for (k=0; k<3; k++) 
	            xc[k] = PeriodicPos((double)Void[iv].Ini[k] + xr[k],LBox[k]); 	  
	  }

	  ic = (int)(xc[0]/GridSize[0]);
	  jc = (int)(xc[1]/GridSize[1]);
	  kc = (int)(xc[2]/GridSize[2]);

	  done = false;
	  p = 0;

	  do {

             MinDist = MaxRadiusSearch/(double)NumShell*(double)(p  );		  
             MaxDist = MaxRadiusSearch/(double)NumShell*(double)(p+1);		  

	     for (in=0; in<NumNeigh[p]; in++) {

	         ii = PeriodicGrid(Neigh[p].i[in] + ic,NumGrid); 
	         jj = PeriodicGrid(Neigh[p].j[in] + jc,NumGrid); 
	         kk = PeriodicGrid(Neigh[p].k[in] + kc,NumGrid); 

		 l = Index1D(ii,jj,kk,NumGrid);

		 if (GridList[l].NumMem == 0) continue;

		 for (m=0; m<GridList[l].NumMem; m++) {
		
	             next = GridList[l].Member[m];		 

		     for (k=0; k<3; k++) 
	                 dx[k] = PeriodicDeltaPos(xc[k] - (double)Tracer[next].Pos[k],LBox[k]);

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
             QSort(SortArr,0,Nsort-1);

	     Radius = 0.5*(SortArr[Nsort-2].val + SortArr[Nsort-1].val);    
   	     Volume = (4.0/3.0)*PI*Radius*Radius*Radius;
	     Delta = (double)(Nsort-1)/Volume/MeanNumTrac - 1.0;

	     if (Delta < DeltaThreshold) {
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
	      Delta = (double)(ir+1)/Volume/MeanNumTrac - 1.0;

              if (Delta < DeltaThreshold && Radius > BiggestRadius) {
		 
		 if (Radius/BiggestRadius - 1.0 >= RadIncrement) CheckRan = 0;
		
	         Void[iv].Rad = (float)Radius;
		 Void[iv].Delta = (float)Delta;
	         Void[iv].Pos[0] = (float)xc[0];	       
	         Void[iv].Pos[1] = (float)xc[1];	       
	         Void[iv].Pos[2] = (float)xc[2];
	         Void[iv].ToF = true;
		 kappa = ir + 1; 
	         BiggestRadius = Radius;
	      
	      } /* Fin lazo Dcum < DeltaThreshold */

	  } /* Fin lazo bines */
	 
	  free(SortArr);
 
      } while (CheckRan < NumRanWalk); /* Fin lazo random */      

      if (Void[iv].ToF) {
	 lambda = (4.0/3.0)*PI*pow((double)Void[iv].Rad,3)*MeanNumTrac;
         Void[iv].Poisson = (double)kappa*log(lambda) - lambda - LnFactorial(kappa); 	 
	 Void[iv].Nran = TotRan;
      }    

  } /* Fin lazo voids */

  for (p=0; p<NumShell; p++) 
      FreeNeighbours(&Neigh[p]);	  
  FreeGridList(GridList,NumGrid);

  StepName.push_back("Finding voids");
  StepTime.push_back(Time(t,OMPcores));

}

void CleanVoids()
{
  int              p,i,j,k,l,m,indx,next,in,it;
  int              ii,jj,kk,ic,jc,kc,NumTrueVoid;
  int              NumNeigh,NumGrid;
  double           dist,Ri,Rj,Vi,Vj,Vij,MinDist,MaxDist;
  double           xi[3],xj[3],dx[3],GridSize[3],GAP;
  struct neighbour Neigh;
  struct sort      *SortArr; 
  struct grid      *GridList;
  clock_t          t;

  fprintf(logfile,"\n CLEANING VOID CATALOGUE BY OVERLAP (TOL = %4.2f) \n",OverlapTol);
  t = clock();

  NumGrid = (int)cbrt((double)NumVoid/10.0);
  GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid)); 
  BuildGridList(GridList,NumGrid,GridSize,1,false);
  
  NumTrueVoid = 0;
  for (i=0; i<NumVoid; i++)
      if (Void[i].ToF) 
	 NumTrueVoid++;	  
  
  fprintf(logfile," | Number of true voids before cleaning = %d (%4.2f %)\n",
		  NumTrueVoid,(double)NumTrueVoid/(double)NumVoid*100.0);
  fflush(logfile);

  SortArr = (struct sort *) malloc(NumTrueVoid*sizeof(struct sort));

  it = 0;
  for (i=0; i<NumVoid; i++) {
      if (Void[i].ToF) {	      
         SortArr[it].val = Void[i].Rad;
         SortArr[it].ord = i;
	 it++;
      }     
  }

  QSort(SortArr,0,NumTrueVoid-1);

  // Selecciono vecinos
  
  GAP = 0.0;
  for (k=0; k<3; k++) 
      if (GridSize[k] > GAP) 
	 GAP = GridSize[k];	  
  GAP *= sqrt(3.0);

  MinDist = 0.0;
  MaxDist = 2.0*MaxRadiusSearch + GAP;  

  SearchNeighbours(&Neigh,&NumNeigh,GridSize,MinDist,MaxDist);
  
  fprintf(logfile," | MinDist - MaxDist = %5.3f - %5.3f [Mpc/h], %d grids \n",MinDist,MaxDist,NumNeigh);
  fflush(logfile);

  for (i=NumTrueVoid-1; i>=0; i--) {

      //if ((NumVoid-i) % 10000 == 0) Progress(NumVoid-i,NumVoid); 

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
    
      for (in=0; in<NumNeigh; in++) {

   	  ii = PeriodicGrid(Neigh.i[in] + ic,NumGrid); 
	  jj = PeriodicGrid(Neigh.j[in] + jc,NumGrid); 
	  kk = PeriodicGrid(Neigh.k[in] + kc,NumGrid); 

	  l = Index1D(ii,jj,kk,NumGrid);

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
	             dx[k] = PeriodicDeltaPos(xi[k] - xj[k],LBox[k]);

                 dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

	         // Volumen de la interseccion, normalizado a Vj
                 // Solo vale para d <= Ri+Rj

	         if (dist > Ri + Rj) continue;

	         Vij  = (PI/12.0/dist);
                 Vij *= pow(Ri + Rj - dist,2);
	         Vij *= (pow(dist,2) + 2.0*dist*(Ri + Rj) - 3.0*pow(Ri - Rj,2));
	         Vij /= Vj;

	         if (Vij > OverlapTol) Void[next].ToF = false;
    
	      } 
    	  }
      }
  }

  NumTrueVoid = 0;
  for (i=0; i<NumVoid; i++)
      if (Void[i].ToF) 
	 NumTrueVoid++;	  
  
  fprintf(logfile," | Number of true void after cleaning = %d (%4.2f %)\n",
		  NumTrueVoid,(double)NumTrueVoid/(double)NumVoid*100.0);

  free(SortArr);
  FreeGridList(GridList,NumGrid);
  FreeNeighbours(&Neigh);

  StepName.push_back("Cleaning voids"); 
  StepTime.push_back(Time(t,1));

}

