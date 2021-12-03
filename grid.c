
#include "allvars.h"
#include "proto.h"

void BuildGridList(struct grid *GridList, int NumGrid, double *GridSize, int who, bool compute_neigh)
{

  int p,i,j,k,l,N,nv,ll;
  int ii,jj,kk,it,jt,kt;

  fprintf(logfile," | Creating grid-list for tracers\n");

  for (k=0; k<3; k++) 
      GridSize[k] = LBox[k]/(double)NumGrid;    

  fprintf(logfile," | Number of grids = %d \n",NumGrid);
  fprintf(logfile," | Grid sizes: x = %f [Mpc/h] \n",GridSize[0]);
  fprintf(logfile," |             y = %f [Mpc/h] \n",GridSize[1]);
  fprintf(logfile," |             z = %f [Mpc/h] \n",GridSize[2]);
  fflush(logfile);

  if (who == 0) 
     N = NumTrac;
  else
     N = NumVoid;

  for (l=0; l<NumGrid*NumGrid*NumGrid; l++)
      GridList[l].NumMem = 0;
   
  for (p=0; p<N; p++) {

      if (who == 0) {	  
         i = (int)(Tracer[p].Pos[0]/GridSize[0]);
         j = (int)(Tracer[p].Pos[1]/GridSize[1]);
         k = (int)(Tracer[p].Pos[2]/GridSize[2]);
      } else {
         i = (int)(Void[p].Pos[0]/GridSize[0]);
         j = (int)(Void[p].Pos[1]/GridSize[1]);
         k = (int)(Void[p].Pos[2]/GridSize[2]);
      }

      if (i == NumGrid) i--;
      if (j == NumGrid) j--;
      if (k == NumGrid) k--;

      l = Index1D(i,j,k,NumGrid);

      GridList[l].NumMem++;
  }

  for (l=0; l<NumGrid*NumGrid*NumGrid; l++) {
      GridList[l].Member = (int *) malloc(GridList[l].NumMem*sizeof(int));	  
      GridList[l].NumMem = 0;
  }

  for (p=0; p<N; p++) {

      if (who == 0) {	  
         i = (int)(Tracer[p].Pos[0]/GridSize[0]);
         j = (int)(Tracer[p].Pos[1]/GridSize[1]);
         k = (int)(Tracer[p].Pos[2]/GridSize[2]);
      } else {
         i = (int)(Void[p].Pos[0]/GridSize[0]);
         j = (int)(Void[p].Pos[1]/GridSize[1]);
         k = (int)(Void[p].Pos[2]/GridSize[2]);
      }

      if (i == NumGrid) i--;
      if (j == NumGrid) j--;
      if (k == NumGrid) k--;

      l = Index1D(i,j,k,NumGrid);

      GridList[l].Member[GridList[l].NumMem] = p;
      GridList[l].NumMem++;

  }

  if (!compute_neigh) return;

  for (l=0; l<NumGrid*NumGrid*NumGrid; l++) 
      GridList[l].Neighbour = (int *) malloc(27*sizeof(int));	  

  for (i=0; i<NumGrid; i++) {
      for (j=0; j<NumGrid; j++) {
          for (k=0; k<NumGrid; k++) {

	      nv = 0;
              l = Index1D(i,j,k,NumGrid);

	      for (it=i-1; it<=i+1; it++) {
		  ii = PeriodicGrid(it,NumGrid);

	          for (jt=j-1; jt<=j+1; jt++) {
		      jj = PeriodicGrid(jt,NumGrid);    

	              for (kt=k-1; kt<=k+1; kt++) {
		          kk = PeriodicGrid(kt,NumGrid);    
              
			  ll = Index1D(ii,jj,kk,NumGrid);
			  GridList[l].Neighbour[nv] = ll;
		          nv++;

		      }
		  }
	      }	  
	  }
      }
  }

  return;
}

void SearchNeighbours(struct neighbour *Neigh, int *NumNeigh, double *GridSize, double MinDist, double MaxDist)
{
  int    i,j,k,NG[3];
  double dist,x,y,z;

  (*NumNeigh) = 0;

  for (k=0; k<3; k++)
      NG[k] = (int)ceil(MaxDist/GridSize[k]);	

  for (i=-NG[0]; i<NG[0]; i++) {
      x = (double)i*GridSize[0];	 

      for (j=-NG[1]; j<NG[1]; j++) {
          y = (double)j*GridSize[1];	  

          for (k=-NG[2]; k<NG[2]; k++) {
              z = (double)k*GridSize[2];	  

              dist = sqrt(x*x + y*y + z*z); 	 

              if (dist >= MinDist && dist <= MaxDist) {
            	 (*Neigh).i.push_back(i);     
           	 (*Neigh).j.push_back(j);     
           	 (*Neigh).k.push_back(k);
                 (*NumNeigh)++;	 
              }
          }
      }
  }  

}

void FreeGridList(struct grid *GridList, int NG)
{
   int i;

   for (i=0; i<NG*NG*NG; i++)  
       free(GridList[i].Member);
   free(GridList);
}

void FreeNeighbours(struct neighbour *Neigh)
{

   (*Neigh).i.clear();
   (*Neigh).j.clear();
   (*Neigh).k.clear();

}

