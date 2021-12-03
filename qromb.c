
#include "allvars.h"
#include "proto.h"

#define EPS   1.0e-6
#define JMAX  20
#define JMAXP (JMAX+1)
#define K 5
#define FUNC(x) ((*func)(x))
#define NR_END 1
#define FREE_ARG char*

double E(double z)
{
  double q;

  q = Cosmo.OmM*(1.0 + z)*(1.0 + z)*(1.0 + z)
    + Cosmo.OmK*(1.0 + z)*(1.0 + z) 
    + Cosmo.OmL;

  q = sqrt(q);

  return(q);
}

double ComovingDistance_integ(double z) 
{
  double q;

  q = 1.0/E(z);

  return(q);
}

double ComovingDistance(double z)
{
  double distance;

  distance = (CVEL/Cosmo.Hub)*qromb(ComovingDistance_integ,0.0,z);

  return(distance);
}

double AngularDistance(double z)
{
  double dcom,dang,dh;

  dcom = ComovingDistance(z);
  dh = CVEL/Cosmo.Hub;

  if (Cosmo.OmK == 0.0) {
     dang = dcom;
  } else if (Cosmo.OmK > 0.0) {
     dang = (dh/sqrt(Cosmo.OmK))*sinh(sqrt(Cosmo.OmK)*dcom/dh);	  
  } else if (Cosmo.OmK < 0.0) {
     dang = (dh/sqrt(fabs(Cosmo.OmK)))*sin(sqrt(fabs(Cosmo.OmK))*dcom/dh);	  
  }

  dang /= (1.0 + z);

  return(dang);  
}	

double *vec(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) {
	   fprintf(stdout,"Numerical Recipes run-time error...\n");
	   fprintf(stdout,"allocation failure in vec()\n");
	   fprintf(stdout,"...now exiting to system...\n");
	   fflush(stdout);
	   exit(EXIT_FAILURE);
	}
	return v-nl+NR_END;
}

void free_vec(double *v, long nl, long nh)
/* free a float vector allocated with vec() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

double trapzd(double (*func)(double), double a, double b, int n)
{
	double x,tnm,sum,del;
	static float s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  double *c,*d;
  //printf("polint ----"); fflush(stdout);

  dif=fabs(x-xa[1]);
  c=vec(1,n);
  d=vec(1,n);
  for (i=1;i<=n;i++) 
  {
    if ( (dift=fabs(x-xa[i])) < dif) 
    {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) 
  {
    for (i=1;i<=n-m;i++) 
    {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ((den=ho-hp) == 0.0) {
	 fprintf(stdout,"Numerical Recipes run-time error...\n");
	 fprintf(stdout,"Error in routine polint\n");
	 fprintf(stdout,"...now exiting to system...\n");
	 fflush(stdout);
	 exit(EXIT_FAILURE);
      }
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free_vec(d,1,n);
  free_vec(c,1,n);
}

double qromb(double (*func)(double), double a, double b)
{
  void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
  double trapzd(double (*func)(double), double a, double b, int n);
  double ss,dss;
  double s[JMAXP+1],h[JMAXP+1];
  int j;

  h[1]=1.0;
  for (j=1;j<=JMAX;j++) 
  {
    s[j]=trapzd(func,a,b,j);
    if (j >= K) 
    {
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) < EPS*fabs(ss)) return ss;
    }
    s[j+1]=s[j];
    h[j+1]=0.25*h[j];
  }

  fprintf(stdout,"Numerical Recipes run-time error...\n");
  fprintf(stdout,"Too many steps in routine qromb\n");
  fprintf(stdout,"...now exiting to system...\n");
  fflush(stdout);
  exit(EXIT_FAILURE);
  return 0.0;
}

void locate(double xx[], int n, double x, int *j)
{
  int ju,jm,jl;
  int ascnd;
  
  jl = 0   ;
  ju = n+1 ;
  ascnd=(xx[n-1] > xx[0]);
  while (ju-jl > 1) 
  {
    jm=(ju+jl) >> 1;
    if (x > xx[jm-1] == ascnd)
      jl = jm ;
    else
      ju = jm ;
  }
  *j = jl ;
}

#undef FUNC
