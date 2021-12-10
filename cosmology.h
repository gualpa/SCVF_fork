
#define WORKSIZE 100000

struct cosmoparam {
  double OmM; 
  double OmL; 
  double OmK;
  double Hub;
};

double E(double, struct cosmoparam *);
double ComovingDistance_integ(double, void *); 
double ComovingDistance(double, struct cosmoparam *);
double AngularDistance(double, struct cosmoparam *);

