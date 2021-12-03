
#include "allvars.h"
#include "proto.h"

int PeriodicGrid(int i, int N)
{
  int ip;
	
  if (i >= N) {
     ip = i - N;	  
  } else if (i < 0) {
     ip = i + N;	  
  } else {
     ip = i;	  
  }

  return ip;
}	

double PeriodicPos(double x, double l)
{
  double xp;
	
  if (x >= l) {
     xp = x - l;	  
  } else if (x < 0.0) {
     xp = x + l;	  
  } else {
     xp = x;	  
  }

  return xp;
}

double PeriodicDeltaPos(double dx, double l)
{
  double dxp;
  double halfbox = 0.5*l;

  if (dx > halfbox) {
     dxp = dx - l;	  
  } else if (dx < -halfbox) {
     dxp = dx + l;	  
  } else {
     dxp = dx;	  
  }

  return dxp;
}

int Index1D(int i, int j, int k, int N)
{
   int indx;
   indx = (k*N + j)*N + i;
   return indx;   
}

vector <int> Index3D(int l, int N)
{
   int          i,j,k;
   vector <int> indx;

   indx.clear();
   i = l % N;
   j = (l - i) / N % N;
   k = ((l - i) / N - j) / N;

   indx.push_back(i);
   indx.push_back(j);
   indx.push_back(k);

   return indx;
}

int CountLines(char *filename)
{
  FILE *fp;
  int  count = 0;
  char string[256];

  fp = SafeOpen(filename, "r");
  while (fgets(string,256,fp)) count++;
  fclose(fp);

  return count; 
}
    
double Time(clock_t ti, int N)
{
  clock_t tf;
  double  tseg,tmin;

  tf = clock();

  tseg = (double)(tf - ti)/(double)CLOCKS_PER_SEC/(double)N;
  tmin = tseg/60.0;

  fprintf(logfile," | -> time taken: %f seg (%f min) \n",tseg,tmin);
 
  fflush(logfile);

  return tseg;
}

void StatsTime()
{
  int     i;	
  double  total;

  fprintf(logfile,"\n TIME STATISTICS\n");

  total = 0.0;
  for (i=0; i<StepTime.size(); i++) 
      total += StepTime[i];
  
  fprintf(logfile," | Total time taken: %f seg (%f min) \n",total,total/60.);
  
  for (i=0; i<StepTime.size(); i++) 
      fprintf(logfile," | %30s: %5.2f %s \n",StepName[i].c_str(),StepTime[i]/total*100.,"%");      
 
  fflush(logfile);

}

void Progress(int par, int tot)
{
  float prog;

  prog = (float)par/(float)tot*100.0;
  fprintf(stdout," | %4.1f %s \r",prog,"%");
  fflush(stdout);

}

double RandomNumber()
{
   double x;
   x = (double)rand()/RAND_MAX;
   return x;
}

double LnFactorial(int n)
{
  int    i;
  double f;

  if (n > 1) {
     f = 0.0;
     for (i=2; i<=n; i++) 
         f += log((double)i);  	  
  } else {
     f = 0.0;	  
  }
	  
  return f;
}

FILE *SafeOpen(char *fname, const char *mode)
{
  FILE *f;

  if (!(f = fopen(fname,mode))) {
     fprintf(stdout,"\nError. File '%s' not found. \n\n",fname);	  
     fflush(stdout);
     exit(EXIT_FAILURE);
  }   

  return f; 
}

