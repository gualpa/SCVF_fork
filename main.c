
#include "allvars.h"
#include "proto.h"

int main(int argc, char **argv) 
{
 
   clock_t t;

   t = clock();   

   if (argc < 2) {
       fprintf(stdout, "\n Error. Missing input file.\n");
       fprintf(stdout, "./main.x <input_param>\n\n");
       exit(EXIT_FAILURE);
   }

   ReadInputFile(argv[1]);

   if (Redshift == 0.0 && GDist == 1) {
      fprintf(stdout,"\nError. Geometrical distortions not available for z = 0\n");
      exit(EXIT_FAILURE);
   }

   omp_set_num_threads(OMPcores);
   fprintf(stdout,"\n ====>>>> Void finder runnning in %d core(s) <<<<==== \n",OMPcores);

   fprintf(stdout,"\nReading tracers... ");fflush(stdout);
   ReadTracers();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing Voronoi tessellation... ");fflush(stdout);
   Voronoi();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nSearching candidates... ");fflush(stdout);
   FindCenters(); 
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nPerforming void identification... ");fflush(stdout);
   FindVoids();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nCleaning void catalogue... ");fflush(stdout);
   CleanVoids();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing void velocities... ");fflush(stdout);
   ComputeVelocity();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing void profiles... ");fflush(stdout);
   ComputeProfiles();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nWrinting void catalogue... ");fflush(stdout);
   WriteVoids();
   fprintf(stdout,"Done.\n\n");fflush(stdout);

   StatsTime();
   
   free(Tracer);
   Void.clear();
   fclose(logfile);

   return(0);
}
