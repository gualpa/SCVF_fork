
#define WORKSIZE 100000
#define GSL_EPS 1.0e-8

struct cosmology {
  double OmM; 
  double OmL; 
  double OmK;
  double Hub;
};

double E(double, struct cosmology *);
double ComovingDistance_integ(double, void *); 
double ComovingDistance(double, struct cosmology *);
double AngularDistance(double, struct cosmology *);

