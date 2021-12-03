
#define EPS   1.0e-6
#define JMAX  20
#define JMAXP (JMAX+1)
#define K 5
#define NR_END 1
#define FREE_ARG char*

double E(double);
double ComovingDistance_integ(double); 
double ComovingDistance(double);
double AngularDistance(double);
double *vec(long nl, long nh);
void   free_vec(double *v, long nl, long nh);
double trapzd(double (*func)(double), double a, double b, int n);
void   polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double qromb(double (*func)(double), double a, double b);
void   locate(double xx[], int n, double x, int *j);

