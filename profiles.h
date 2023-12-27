
struct profile {
   float DeltaDiff;
   float DeltaCum;
   float Velocity;
   float Ri,Rm,Rs;
};
varConfiguration compute_profiles(varConfiguration VarConfigAux);

//void compute_profiles();
//void bin2ascii_profile(int voidID);
varConfiguration bin2ascii_profile(int voidID, varConfiguration VarConfigAux);
