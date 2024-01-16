/**
 * @file allvars.h
 * @brief Programa que identifica voids astronómicos
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
#include "io.h"
#include "voronoi.h" 
#include "finder.h"
#include "velocity.h"
#include "profiles.h"
#include "tools.h"

/**
 * @brief Punto de entrada principal del programa astronómico para identificar y analizar voids.
 *
 * Esta función principal ejecuta un programa para identificar voids astronómicos basado en ciertos parámetros
 * de entrada. Lee un archivo de parámetros de entrada, procesa los trazadores astronómicos, realiza cálculos
 * de Voronoi, encuentra candidatos a voids, identifica voids, limpia el catálogo de voids, calcula velocidades
 * de voids, calcula perfiles de voids y escribe un catálogo de voids en un archivo.
 *
 * @param argc Cantidad de argumentos de línea de comandos.
 * @param argv Nombre del archivo de configuración con los parámetros astronómicos. Ejemplo "VF.param".
 * @return Entero que indica el estado de salida del programa (0 para éxito, otro valor para error).
 */
int main(int argc, char **argv) 
{
   varConfiguration VarConfigAux;
   logs LogAux;
   vector <tracers> TracerAux;
   vector <voids> VoidAux;

   if (argc < 2) {
       fprintf(stdout, "\n Error. Missing input file and flags.\n");
       fprintf(stdout, "./main.x <input_param> [<run flag>] [...] \n\n");
       exit(EXIT_FAILURE);
   } else if (argc == 2) {
       VarConfigAux.RunFlag = 0;
   } else {
       sscanf(argv[2], "%d", &VarConfigAux.RunFlag);
   }

   VarConfigAux = read_input_file(argv[1],VarConfigAux,LogAux);

   if (VarConfigAux.RunFlag == 1) {
      if (argc < 4) {
	 fprintf(stdout, "\n Error. Missing void ID.\n");
         fprintf(stdout, "./main.x <input_param> 1 <void_ID> \n\n");
         exit(EXIT_FAILURE);
      }	      
      int voidID; 
      sscanf(argv[3], "%d", &voidID);
      VarConfigAux = bin2ascii_profile(voidID, VarConfigAux, LogAux);
      exit(EXIT_SUCCESS);
   } 

   if (VarConfigAux.Redshift == 0.0 && VarConfigAux.GDist == 1) {
      fprintf(stdout,"\nError. Geometrical distortions not available for z = 0\n");
      exit(EXIT_FAILURE);
   }

   omp_set_num_threads(VarConfigAux.OMPcores);
   fprintf(stdout,"\n ====>>>> Void finder runnning in %d core(s) <<<<==== \n",VarConfigAux.OMPcores);

   fprintf(stdout,"\nReading tracers... ");fflush(stdout);
   TracerAux = read_tracers(VarConfigAux, LogAux);
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing Voronoi tessellation... ");fflush(stdout);
   compute_voronoi(VarConfigAux, LogAux, TracerAux);
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nSearching candidates... ");fflush(stdout);
   VoidAux = find_void_candidates(VarConfigAux, LogAux, TracerAux);
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nPerforming void identification... ");fflush(stdout);
   VoidAux = find_voids(VarConfigAux, LogAux, TracerAux, VoidAux);
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nCleaning void catalogue... ");fflush(stdout);
   VoidAux = clean_voids(VarConfigAux, LogAux, VoidAux);
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing void velocities... ");fflush(stdout);
   compute_velocity(VarConfigAux, LogAux, TracerAux, VoidAux);
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing void profiles... ");fflush(stdout);
   compute_profiles(VarConfigAux, LogAux, TracerAux, VoidAux);
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nWrinting void catalogue... ");fflush(stdout);
   write_voids(VarConfigAux, LogAux, VoidAux);
   fprintf(stdout,"Done.\n\n");fflush(stdout);

   time_resume(VarConfigAux, LogAux);
   
   TracerAux.clear();
   VoidAux.clear();
   fclose(LogAux.logfile);

   return(0);
}
