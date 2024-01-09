
struct profile {
   float DeltaDiff;
   float DeltaCum;
   float Velocity;
   float Ri,Rm,Rs;
};
void compute_profiles(varConfiguration &VarConfigAux, logs &LogAux, vector <tracers> TracerAux, vector <voids> &VoidAux);
varConfiguration bin2ascii_profile(int voidID, varConfiguration VarConfigAux, logs &LogAux);
