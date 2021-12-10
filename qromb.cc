
#include "allvars.h"
#include "qromb.h"
#include "gsl_integration.h"
#include "gsl_math.h"

double E(double z, struct cosmology *C)
{
  double q;

  q = C->OmM*(1.0 + z)*(1.0 + z)*(1.0 + z)
    + C->OmK*(1.0 + z)*(1.0 + z) 
    + C->OmL;

  q = sqrt(q);

  return(q);
}

double ComovingDistance_integ(double z, void *p) 
{
  double q;
  struct cosmology *C = (struct cosmology *)p;

  q = 1.0/E(z,C);
  
  return(q);
}

double ComovingDistance(double z, struct cosmology *C)
{
  double distance,result,abserr;
  gsl_function F;
  gsl_integration_workspace *workspace;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &ComovingDistance_integ;
  F.params = C;
  
  gsl_integration_qag(&F,0.0,z,0,1.0e-8,WORKSIZE,GSL_INTEG_GAUSS41,workspace,&result,&abserr);

  distance = (CVEL/C->Hub)*result;

  gsl_integration_workspace_free(workspace);

  return(distance);
}

double AngularDistance(double z, struct cosmology *C)
{
  double dcom,dang,dh;

  dcom = ComovingDistance(z,C);
  dh = CVEL/C->Hub;

  if (C->OmK == 0.0) {
     dang = dcom;
  } else if (C->OmK > 0.0) {
     dang = (dh/sqrt(C->OmK))*sinh(sqrt(C->OmK)*dcom/dh);	  
  } else if (C->OmK < 0.0) {
     dang = (dh/sqrt(fabs(C->OmK)))*sin(sqrt(fabs(C->OmK))*dcom/dh);	  
  }

  dang /= (1.0 + z);

  return(dang);  
}	

