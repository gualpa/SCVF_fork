#include "allvars.h"
varConfiguration read_input_file(char *filename, varConfiguration VarConfigAux, logs &LogAux);
void read_tracers(varConfiguration &VarConfigAux, logs &LogAux);
void read_tracers_ascii(varConfiguration &VarConfigAux);
void read_tracers_binary(varConfiguration &VarConfigAux);
void read_tracers_gadget(varConfiguration &VarConfigAux);
void read_tracers_mxxl(varConfiguration &VarConfigAux, logs &LogAux);
void redshift_space_distortions(varConfiguration VarConfigAux, logs &LogAux);
void geometrical_distortions(varConfiguration &VarConfigAux, logs &LogAux);
void write_voids(varConfiguration VarConfigAux, logs &LogAux);