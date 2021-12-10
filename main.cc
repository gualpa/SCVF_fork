
#include "allvars.h"
#include "io.h"
#include "voronoi.h" 
#include "finder.h"
#include "velocity.h"
#include "profiles.h"
#include "tools.h"

int main(int argc, char **argv) 
{
 
   clock_t t;

   t = clock();   

   if (argc < 2) {
       fprintf(stdout, "\n Error. Missing input file.\n");
       fprintf(stdout, "./main.x <input_param>\n\n");
       exit(EXIT_FAILURE);
   }

   read_input_file(argv[1]);

   if (Redshift == 0.0 && GDist == 1) {
      fprintf(stdout,"\nError. Geometrical distortions not available for z = 0\n");
      exit(EXIT_FAILURE);
   }

   omp_set_num_threads(OMPcores);
   fprintf(stdout,"\n ====>>>> Void finder runnning in %d core(s) <<<<==== \n",OMPcores);

   fprintf(stdout,"\nReading tracers... ");fflush(stdout);
   read_tracers();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing Voronoi tessellation... ");fflush(stdout);
   compute_voronoi();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nSearching candidates... ");fflush(stdout);
   find_void_candidates(); 
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nPerforming void identification... ");fflush(stdout);
   find_voids();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nCleaning void catalogue... ");fflush(stdout);
   clean_voids();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing void velocities... ");fflush(stdout);
   compute_velocity();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing void profiles... ");fflush(stdout);
   compute_profiles();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nWrinting void catalogue... ");fflush(stdout);
   write_voids();
   fprintf(stdout,"Done.\n\n");fflush(stdout);

   time_resume();
   
   free(Tracer);
   Void.clear();
   fclose(logfile);

   return(0);
}
