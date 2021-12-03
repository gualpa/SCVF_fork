
/* I/O */

void ReadInputFile(char *);
void ReadTracers(void);
void ReadTracers_ASCII(void);
void ReadTracers_BINARY(void);
void ReadTracers_MXXL(void);
void ReadTracers_GADGET1(void);
void ReadTracers_GADGET2(void);
void WriteVoids(void);

/* Voronoi tessellation */

void Voronoi(void);
void ComputeVoronoi(void);
void ReadVoronoi_ASCII(void);
void WriteVoronoi_ASCII(void);
void ReadVoronoi_MXXL(void);
void WriteVoronoi_MXXL(void);

/* Tools */
void         SearchNeighbours(struct neighbour *, int *, double *, double, double);
void         FreeNeighbours(struct neighbour *);
void         BuildGridList(struct grid *, int, double *, int, bool);
void         FreeGridList(struct grid *, int);
int          Index1D(int,int,int,int);
vector <int> Index3D(int,int);
double       PeriodicPos(double,double);
double       PeriodicDeltaPos(double,double);
int          PeriodicGrid(int,int);
double       Time(clock_t, int);
void         StatsTime();
void         Progress(int,int);
int          CountLines(char *);
double       LnFactorial(int);
double       RandomNumber(void);
FILE*        SafeOpen(char *, const char *); 

/* Sorting */

void QSort(struct sort *, int, int);

/* Finder */

void FindCenters(void);
void FindVoids(void);
void CleanVoids(void);
void ComputeVelocity(void);
void ComputeProfiles(void);

/* Cosmology */

double E(double);
double ComovingDistance_integ(double); 
double ComovingDistance(double);
double AngularDistance(double);
void   GeometricalDistortions(void);
void   RedshiftSpaceDistortions(void);

/* NR */

double *vec(long nl, long nh);
void   free_vec(double *v, long nl, long nh);
double trapzd(double (*func)(double), double a, double b, int n);
void   polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double qromb(double (*func)(double), double a, double b);
void   locate(double xx[], int n, double x, int *j);

