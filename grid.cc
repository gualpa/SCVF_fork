
#include "allvars.h"
#include "tools.h"
#include "grid.h"

template <class T>
void build_grid_list(T Points, int NumPoints, struct grid *GridList, int NumGrid, double *GridSize, bool compute_neigh)
{

  fprintf(logfile," | Creating grid-list \n");

  for (int k=0; k<3; k++) 
      GridSize[k] = LBox[k]/(double)NumGrid;    

  fprintf(logfile," | Number of grids = %d \n",NumGrid);
  fprintf(logfile," | Grid sizes: x = %f [Mpc/h] \n",GridSize[0]);
  fprintf(logfile," |             y = %f [Mpc/h] \n",GridSize[1]);
  fprintf(logfile," |             z = %f [Mpc/h] \n",GridSize[2]);
  fflush(logfile);

  for (int l=0; l<NumGrid*NumGrid*NumGrid; l++)
      GridList[l].NumMem = 0;
   
  for (int p=0; p<NumPoints; p++) {

      int i = (int)(Points[p].Pos[0]/GridSize[0]);
      int j = (int)(Points[p].Pos[1]/GridSize[1]);
      int k = (int)(Points[p].Pos[2]/GridSize[2]);

      if (i == NumGrid) i--;
      if (j == NumGrid) j--;
      if (k == NumGrid) k--;

      int l = index_1d(i,j,k,NumGrid);

      GridList[l].NumMem++;
  }

  for (int l=0; l<NumGrid*NumGrid*NumGrid; l++) {
      GridList[l].Member = (int *) malloc(GridList[l].NumMem*sizeof(int));	  
      GridList[l].NumMem = 0;
  }

  for (int p=0; p<NumPoints; p++) {

      int i = (int)(Points[p].Pos[0]/GridSize[0]);
      int j = (int)(Points[p].Pos[1]/GridSize[1]);
      int k = (int)(Points[p].Pos[2]/GridSize[2]);

      if (i == NumGrid) i--;
      if (j == NumGrid) j--;
      if (k == NumGrid) k--;

      int l = index_1d(i,j,k,NumGrid);

      GridList[l].Member[GridList[l].NumMem] = p;
      GridList[l].NumMem++;

  }

  if (!compute_neigh) return;

  int nv;

  for (int l=0; l<NumGrid*NumGrid*NumGrid; l++) 
      GridList[l].Neighbour = (int *) malloc(27*sizeof(int));	  

  for (int i=0; i<NumGrid; i++) {
      for (int j=0; j<NumGrid; j++) {
          for (int k=0; k<NumGrid; k++) {

	      nv = 0;
              int l = index_1d(i,j,k,NumGrid);

	      for (int it=i-1; it<=i+1; it++) {
		  int ii = periodic_grid(it,NumGrid);

	          for (int jt=j-1; jt<=j+1; jt++) {
		      int jj = periodic_grid(jt,NumGrid);    

	              for (int kt=k-1; kt<=k+1; kt++) {
		          int kk = periodic_grid(kt,NumGrid);    
              
			  int ll = index_1d(ii,jj,kk,NumGrid);
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

void search_neighbours(struct neighbour *Neigh, int *NumNeigh, double *GridSize, double MinDist, double MaxDist)
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

void free_grid_list(struct grid *GridList, int NG)
{
   int i;

   for (i=0; i<NG*NG*NG; i++)  
       free(GridList[i].Member);
   free(GridList);
}

void free_neighbours(struct neighbour *Neigh)
{

   (*Neigh).i.clear();
   (*Neigh).j.clear();
   (*Neigh).k.clear();

}

template void build_grid_list<struct tracers*>(struct tracers* Points, int NumPoints, struct grid *GridList, 
		                             int NumGrid, double *GridSize, bool compute_neigh);
template void build_grid_list<vector <voids>>(vector <voids> Points, int NumPoints, struct grid *GridList, 
		                            int NumGrid, double *GridSize, bool compute_neigh);
