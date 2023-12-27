
int          index_1d(int,int,int,int);
vector <int> index_3d(int,int);
double       periodic_position(double,double);
double       periodic_delta(double,double);
int          periodic_grid(int,int);
double       get_time(clock_t, int, varConfiguration VarConfigAux);
void         time_resume(varConfiguration VarConfigAux);
int          count_lines(char *);
double       ln_factorial(int);
double       random_number(void);
FILE*        safe_open(char *, const char *);
