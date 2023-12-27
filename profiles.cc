
#include "allvars.h"
#include "grid.h"
#include "tools.h"
#include "profiles.h"

varConfiguration compute_profiles(varConfiguration VarConfigAux)
{
   int            i,k,ic,jc,kc,l,ii,jj,kk,next,ibin,in,m,NumGrid;
   double         xc[3],xt[3],dx[3],vt[3],dist,GridSize[3];
   double         dR,MinDist,MaxDist,GAP;
   double         VRad,Vol,DeltaMax;
   float          Radius;
   vector <int>   Indx;
   struct profile Prof[VarConfigAux.NumProfileBins];
   char           TxtFile[MAXCHAR],BinFile[MAXCHAR];
   int            NumQuery;
   struct query   Query;
   struct grid    *GridList;
   FILE           *ftxt,*fbin;
   clock_t        t;

   fprintf(VarConfigAux.logfile,"\n COMPUTING VOID PROFILES \n");
   t = clock();

   NumGrid = (int)(VarConfigAux.BoxSize/VarConfigAux.ProxyGridSize);
   GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid));
   build_grid_list(Tracer,VarConfigAux.NumTrac,GridList,NumGrid,GridSize,false,VarConfigAux);

   // Only for true voids

   for (i=0; i<VarConfigAux.NumVoid; i++)
       if (Void[i].ToF) 
	  Indx.push_back(i);       

   dR = (log10(VarConfigAux.MaxProfileDist)-log10(VarConfigAux.MinProfileDist))/(double)VarConfigAux.NumProfileBins;
   
   // Selecciono grides

   MinDist = 0.0;
   MaxDist = 0.0;
   for (i=0; i<VarConfigAux.NumVoid; i++) {
       if (!Void[i].ToF) continue;
       if (Void[i].Rad > MaxDist) MaxDist = Void[i].Rad;
   }
   MaxDist *= VarConfigAux.MaxProfileDist;
   query_grid(&Query,GridSize,MinDist,MaxDist);
   NumQuery = Query.i.size();
   GAP = 0.5*sqrt(3.0)*max_grid_size(GridSize);
   
   fprintf(VarConfigAux.logfile," | MinDist - MaxDist = %5.3f - %5.3f [Mpc/h], %d grids \n",MinDist,MaxDist,NumQuery);
   fflush(VarConfigAux.logfile);

   if (VarConfigAux.WriteProfiles == 2) {
      sprintf(BinFile,"%s/profiles.bin",VarConfigAux.PathProfiles);
      fbin = safe_open(BinFile,"w");
   }

   #pragma omp parallel for default(none) schedule(dynamic)                   \
    shared(VarConfigAux,Void,Tracer,NumQuery,Query,dR,NumGrid,GridSize,GridList,   \
           Indx,GAP,fbin,BinFile)                       \
    private(i,m,k,ii,jj,kk,l,Radius,ic,jc,kc,xc,xt,dx,vt,next,Prof,dist,VRad, \
	    ibin,DeltaMax,Vol,ftxt,in,TxtFile)

   for (i=0; i<(int)Indx.size(); i++) {
       
       for (k=0; k<VarConfigAux.NumProfileBins; k++) {
           Prof[k].DeltaDiff = 0.0;
           Prof[k].DeltaCum = 0.0;
           Prof[k].Velocity = 0.0;
       }
       
       Radius = Void[Indx[i]].Rad;
       for (k=0; k<3; k++) 
	   xc[k] = Void[Indx[i]].Pos[k];
       Void[Indx[i]].Dtype = 0.0;
      
       ic = (int)(xc[0]/GridSize[0]);
       jc = (int)(xc[1]/GridSize[1]);
       kc = (int)(xc[2]/GridSize[2]);

       for (in=0; in<NumQuery; in++) {

       	   ii = Query.i[in]; 
	   jj = Query.j[in]; 
	   kk = Query.k[in]; 

	   dist = (double)(ii*ii)*(GridSize[0]*GridSize[0])
	        + (double)(jj*jj)*(GridSize[1]*GridSize[1])
	        + (double)(kk*kk)*(GridSize[2]*GridSize[2]);
	   dist = sqrt(dist);

	   if (dist > VarConfigAux.MaxProfileDist*Radius+GAP) continue;

       	   ii = periodic_grid(ii + ic, NumGrid); 
	   jj = periodic_grid(jj + jc, NumGrid); 
	   kk = periodic_grid(kk + kc, NumGrid); 	       
	   
	   l = index_1d(ii,jj,kk,NumGrid);

           if (GridList[l].NumMem == 0) continue;

           for (m=0; m<GridList[l].NumMem; m++) {
		
	       next = GridList[l].Member[m];	

	       for (k=0; k<3; k++) {
	           xt[k] = (double)Tracer[next].Pos[k];	 
	           vt[k] = (double)(Tracer[next].Vel[k] - Void[Indx[i]].Vel[k]);
	           dx[k] = periodic_delta(xt[k] - xc[k],VarConfigAux.LBox[k])/Radius;
	       }

               dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

	       if (dist > VarConfigAux.MinProfileDist && dist < VarConfigAux.MaxProfileDist) {
           
	          ibin = (int)((log10(dist)-log10(VarConfigAux.MinProfileDist))/dR);

	          VRad = vt[0]*dx[0] + vt[1]*dx[1] + vt[2]*dx[2];
	          VRad /= dist;

	          Prof[ibin].Velocity += VRad;
	          Prof[ibin].DeltaDiff += 1.0; 

	       }
	   }
       }
       
       for (k=0; k<VarConfigAux.NumProfileBins; k++) {
	   if (Prof[k].DeltaDiff < 3.0) {
	       Prof[k].DeltaDiff = 0.0;
               Prof[k].Velocity = 0.0;
           } else {
               Prof[k].Velocity /= Prof[k].DeltaDiff;	   
	   }
           for (kk=0; kk<=k; kk++) 
	       Prof[k].DeltaCum += Prof[kk].DeltaDiff;
       }

       DeltaMax = -1.0;
       for (k=0; k<VarConfigAux.NumProfileBins; k++) {

	   Prof[k].Ri = (float)(k    )*dR + log10(VarConfigAux.MinProfileDist);
	   Prof[k].Rm = (float)(k+0.5)*dR + log10(VarConfigAux.MinProfileDist);
	   Prof[k].Rs = (float)(k+1.0)*dR + log10(VarConfigAux.MinProfileDist);

	   Prof[k].Ri = pow(10.0,Prof[k].Ri)*Radius;
	   Prof[k].Rm = pow(10.0,Prof[k].Rm)*Radius;
	   Prof[k].Rs = pow(10.0,Prof[k].Rs)*Radius;

	   Vol = (4.0/3.0)*PI*(pow(Prof[k].Rs,3) - pow(Prof[k].Ri,3));
	   Prof[k].DeltaDiff = Prof[k].DeltaDiff/Vol/VarConfigAux.MeanNumTrac - 1.0;

	   Vol = (4.0/3.0)*PI*pow(Prof[k].Rs,3);
	   Prof[k].DeltaCum = Prof[k].DeltaCum/Vol/VarConfigAux.MeanNumTrac - 1.0;

	   if (Prof[k].Rs < 2.0*Radius || Prof[k].Rs > 3.0*Radius) continue;

	   if (Prof[k].DeltaCum > DeltaMax) DeltaMax = Prof[k].DeltaCum;

       }

       Void[Indx[i]].Dtype = DeltaMax;
       
       if (VarConfigAux.WriteProfiles == 1) {
          sprintf(TxtFile,"%s/profile_void_%d.dat",VarConfigAux.PathProfiles,i);
          ftxt = safe_open(TxtFile,"w");
          for (k=0; k<VarConfigAux.NumProfileBins; k++)
	      fprintf(ftxt,"%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n",
			   Prof[k].Ri,Prof[k].Rm,Prof[k].Rs,Prof[k].DeltaDiff,
			   Prof[k].DeltaCum,Prof[k].Velocity,Radius);  
          fclose(ftxt);
       }
       
       if (VarConfigAux.WriteProfiles == 2) {
          #pragma omp critical
	  {
	     fwrite(&i,sizeof(int),1,fbin);
             fwrite(&Radius,sizeof(float),1,fbin);
             fwrite(&VarConfigAux.NumProfileBins,sizeof(int),1,fbin);
	     for (k=0; k<VarConfigAux.NumProfileBins; k++) {
	         fwrite(&Prof[k].Ri,sizeof(float),1,fbin);
	         fwrite(&Prof[k].Rm,sizeof(float),1,fbin);
	         fwrite(&Prof[k].Rs,sizeof(float),1,fbin);
	         fwrite(&Prof[k].DeltaDiff,sizeof(float),1,fbin);
	         fwrite(&Prof[k].DeltaCum,sizeof(float),1,fbin);
	         fwrite(&Prof[k].Velocity,sizeof(float),1,fbin);
	     }
	  }
       }

       Void[Indx[i]].Dtype = DeltaMax;
   } 
   
   if (VarConfigAux.WriteProfiles == 2) fclose(fbin);

   Indx.clear();
   free_query_grid(&Query);
   free_grid_list(GridList,NumGrid);
   
   VarConfigAux.StepName.push_back("Computing profiles");
   VarConfigAux.StepTime.push_back(get_time(t,VarConfigAux.OMPcores,VarConfigAux));
   return VarConfigAux;
}

varConfiguration bin2ascii_profile(int voidID, varConfiguration VarConfigAux)
{
   int            SkipBlock,iv,k;
   float          Radius;  
   struct profile Prof[VarConfigAux.NumProfileBins];
   FILE           *fbin,*fout;
   char           filename[MAXCHAR];
   
   SkipBlock = sizeof(float) + sizeof(int) + 6*VarConfigAux.NumProfileBins*sizeof(float);
   sprintf(filename,"%s/profiles.bin",VarConfigAux.PathProfiles);
   fbin = safe_open(filename,"r");

   do {
   
      fread(&iv,sizeof(int),1,fbin);

      if (iv+1 != voidID)  
         fseek(fbin,SkipBlock,SEEK_CUR);
      else {
         fread(&Radius,sizeof(float),1,fbin);
         fread(&VarConfigAux.NumProfileBins,sizeof(int),1,fbin);
	 for (k=0; k<VarConfigAux.NumProfileBins; k++) {
	     fread(&Prof[k].Ri,sizeof(float),1,fbin);
	     fread(&Prof[k].Rm,sizeof(float),1,fbin);
	     fread(&Prof[k].Rs,sizeof(float),1,fbin);
	     fread(&Prof[k].DeltaDiff,sizeof(float),1,fbin);
	     fread(&Prof[k].DeltaCum,sizeof(float),1,fbin);
	     fread(&Prof[k].Velocity,sizeof(float),1,fbin);
	 }
      }
   } while (iv+1 != voidID);
   fclose(fbin);

   fprintf(VarConfigAux.logfile,"\n Writting profile for void %d, with radius %5.3f [Mpc/h] \n",voidID,Radius);

   sprintf(filename,"%s/profile_void_%d.dat",VarConfigAux.PathProfiles,voidID);
   fout = safe_open(filename,"w");

   for (k=0; k<VarConfigAux.NumProfileBins; k++)
       fprintf(fout,"%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n",
		     Prof[k].Ri,Prof[k].Rm,Prof[k].Rs,Prof[k].DeltaDiff,
		     Prof[k].DeltaCum,Prof[k].Velocity,Radius);  
   fclose(fout);
   return VarConfigAux;
}
