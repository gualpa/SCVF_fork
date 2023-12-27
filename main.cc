
#include "allvars.h"
#include "io.h"
#include "voronoi.h" 
#include "finder.h"
#include "velocity.h"
#include "profiles.h"
#include "tools.h"

int main(int argc, char **argv) 
{
   varConfiguration VarConfigAux;
   if (argc < 2) {
       fprintf(stdout, "\n Error. Missing input file and flags.\n");
       fprintf(stdout, "./main.x <input_param> [<run flag>] [...] \n\n");
       exit(EXIT_FAILURE);
   } else if (argc == 2) {
       VarConfigAux.RunFlag = 0;
   } else {
       sscanf(argv[2], "%d", &VarConfigAux.RunFlag);
   }

   VarConfigAux = read_input_file(argv[1],VarConfigAux);

   if (VarConfigAux.RunFlag == 1) {
      if (argc < 4) {
	 fprintf(stdout, "\n Error. Missing void ID.\n");
         fprintf(stdout, "./main.x <input_param> 1 <void_ID> \n\n");
         exit(EXIT_FAILURE);
      }	      
      int voidID; 
      sscanf(argv[3], "%d", &voidID);
      VarConfigAux = bin2ascii_profile(voidID, VarConfigAux);
      exit(EXIT_SUCCESS);
   } 

   if (VarConfigAux.Redshift == 0.0 && VarConfigAux.GDist == 1) {
      fprintf(stdout,"\nError. Geometrical distortions not available for z = 0\n");
      exit(EXIT_FAILURE);
   }

   omp_set_num_threads(VarConfigAux.OMPcores);
   fprintf(stdout,"\n ====>>>> Void finder runnning in %d core(s) <<<<==== \n",VarConfigAux.OMPcores);

   fprintf(stdout,"\nReading tracers... ");fflush(stdout);
   VarConfigAux = read_tracers(VarConfigAux);
   fprintf(stdout,"Done.\n");fflush(stdout);

 fprintf(stdout,"3  VarConfigAux.NumTrac %i  \n",VarConfigAux.NumTrac);    fflush(stdout);
   fprintf(stdout,"\nComputing Voronoi tessellation... ");fflush(stdout);
   VarConfigAux = compute_voronoi(VarConfigAux);
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nSearching candidates... ");fflush(stdout);
   VarConfigAux = find_void_candidates(VarConfigAux);
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nPerforming void identification... ");fflush(stdout);
   VarConfigAux = find_voids(VarConfigAux);
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nCleaning void catalogue... ");fflush(stdout);
   VarConfigAux = clean_voids(VarConfigAux);
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing void velocities... ");fflush(stdout);
   VarConfigAux = compute_velocity(VarConfigAux);
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing void profiles... ");fflush(stdout);
   VarConfigAux = compute_profiles(VarConfigAux);
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nWrinting void catalogue... ");fflush(stdout);
   write_voids(VarConfigAux);
   fprintf(stdout,"Done.\n\n");fflush(stdout);

   time_resume(VarConfigAux);
   
   Tracer.clear();
   Void.clear();
   fclose(VarConfigAux.logfile);

   return(0);
}
