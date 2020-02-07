/**************************************************************************/
/**************************************************************************/
/* 
   This is the Runge-Kutta driver program using adaptive stepsize control. 
   It integrates, starting with the values ystart[1..N], from x1 to x2 with accuracy eps. 
   Intermediate values are stored in the global arrays. 
   The initial stepsize guess is h1, with hmin the minimum stepsize (can be zero). 
   Upon completion, nok and nbad store the number of good and bad (but retried and  xed) steps. 
   The  nal values are stored in ystart. We assume a routine called derivs which calculates 
   the derivatives of yp and a routine called rkqs which performs a single Runge-Kutta integration step.
*/
/**************************************************************************/
/**************************************************************************/
// 6/8/2003   This routine has been modified by A. Bueno, since it was written 
//            in a highly unacceptable FORTRAN influenced style, taken as it is 
//            from the Numerical Recipes

#include "odeint.h"
#define NRANSI 
#include "nrutil.h" 


const int kmax=MAXSTEPS; /* defined in odeint.c */
//extern int kmax;
extern int kount; 
extern double *xp,**yp,dxsav; 

void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1, double hmin, int *nok, int *nbad, void (*derivs)(double, double [], double []), void (*rkqs)(double [], double [], int, double *, double, double, double [], double *, double *, void (*)(double, double [], double []))) 
{ 
  int nstp,i; 
  double xsav,x,hnext,hdid,h; 
  double *yscal,*y,*dydx; 

  int size=nvar*sizeof(double);
  yscal=(double *)malloc(size); 
  y=(double *)malloc(size); 
  dydx=(double *)malloc(size); 

  //  yscal=vector(1,nvar); 
  //  y=vector(1,nvar); 
  //  dydx=vector(1,nvar); 


  x=x1; h=SIGN(h1,x2-x1); 
  *nok = (*nbad) = kount = 0; 
  for (i=0;i<nvar;i++) 
    y[i]=ystart[i]; 
  if (kmax > 0) xsav=x-dxsav*2.0; 
  for (nstp=0;nstp<MAXSTP;nstp++) 
    { (*derivs)(x,y,dydx); 
    for (i=0;i<nvar;i++) 
      yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY; 
    if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) 
      { xp[kount]=x;  // remove ++ from kount 
      for (i=0;i<nvar;i++) 
	yp[i][kount++]=y[i]; //add ++ to kount
      xsav=x; 
      } 
    if ((x+h-x2)*(x+h-x1) > 0.0) 
      h=x2-x; 
    (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs); 
    if (hdid == h) 
      ++(*nok); 
    else 
      ++(*nbad); 
    if ((x-x2)*(x2-x1) >= 0.0) 
      { for (i=0;i<nvar;i++) 
	ystart[i]=y[i]; 
      if (kmax) { 
	xp[kount]=x; // remove ++ from kount 
	for (i=0;i<nvar;i++) 
	  yp[i][kount++]=y[i]; //add ++ to kount
      } 
      free(dydx); 
      free(y); 
      free(yscal); 
      //      free_vector(dydx,1,nvar); 
      //      free_vector(y,1,nvar); 
      //      free_vector(yscal,1,nvar); 
      return; 
     } 
    if (fabs(hnext) <= hmin) 
      nrerror((char*)("Step size too small in odeint")); 
    h=hnext; 
    } 
  nrerror((char*)("Too many steps in routine odeint")); 
} 
#undef MAXSTP 
#undef TINY 
//#undef NRANSI

/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/
#define SAFETY 0.9 
#define PGROW -0.2 
#define PSHRNK -0.25 
#define ERRCON 1.89e-4 

void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps, double yscal[], double *hdid, double *hnext, void (*derivs)(double, double [], double [])) 
{ 
  void rkck(double y[], double dydx[], int n, double x, double h, double yout[], double yerr[], void (*derivs)(double, double [], double [])); 
  int i; 
  double errmax,h,xnew,*yerr,*ytemp; 

  int size=n*sizeof(double);
  yerr=(double *)malloc(size); 
  ytemp=(double *)malloc(size); 
  //  yerr=vector(1,n); 
  //  ytemp=vector(1,n); 

  h=htry; 
  for (;;) { 
    rkck(y,dydx,n,*x,h,ytemp,yerr,derivs); 
    errmax=0.0; 
    for (i=0;i<n;i++) 
      errmax=FMAX(errmax,fabs(yerr[i]/yscal[i])); 
    errmax /= eps; 
    if (errmax > 1.0) { 
      h=SAFETY*h*pow(errmax,PSHRNK); 
      if (h < 0.1*h) 
	h *= 0.1; 
      xnew=(*x)+h; 
      if (xnew == *x) 
	nrerror((char*)("stepsize underflow in rkqs")); 
      continue; 
    } 
    else { 
      if (errmax > ERRCON) 
	*hnext=SAFETY*h*pow(errmax,PGROW);
      else *hnext=5.0*h; 
      *x += (*hdid=h); 
      for (i=0;i<n;i++) 
	y[i]=ytemp[i]; 
      break; 
    } 
  } 

  free(ytemp);
  free(yerr);
  //  free_vector(ytemp,1,n); 
  //  free_vector(yerr,1,n); 
} 
#undef SAFETY 
#undef PGROW 
#undef PSHRNK 
#undef ERRCON 
//#undef NRANSI
/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/
void rkck(double y[], double dydx[], int n, double x, double h, double yout[], double yerr[], void (*derivs)(double, double [], double [])) 
{ 
  int i; 
  static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2, b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2, b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0, b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0, b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0, c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0, dc5 = -277.0/14336.0; 
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0, dc4=c4-13525.0/55296.0,dc6=c6-0.25; 
  double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp; 

  int size=n*sizeof(double);

  ak2=(double*) malloc(size);
  ak3=(double*) malloc(size);
  ak4=(double*) malloc(size);
  ak5=(double*) malloc(size);
  ak6=(double*) malloc(size);
  ytemp=(double*) malloc(size);
  //  ak2=vector(1,n); 
  //  ak3=vector(1,n); 
  //  ak4=vector(1,n); 
  //  ak5=vector(1,n); 
  //  ak6=vector(1,n); 
  //  ytemp=vector(1,n); 
  for (i=0;i<n;i++) 
    ytemp[i]=y[i]+b21*h*dydx[i]; 
  (*derivs)(x+a2*h,ytemp,ak2); 
  for (i=0;i<n;i++) 
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]); 
  (*derivs)(x+a3*h,ytemp,ak3); 
  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]); 
  (*derivs)(x+a4*h,ytemp,ak4); 
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]); 
  (*derivs)(x+a5*h,ytemp,ak5); 
  for (i=0;i<n;i++) 
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]); 
  (*derivs)(x+a6*h,ytemp,ak6); 
  for (i=0;i<n;i++) 
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]); 
  for (i=0;i<n;i++) 
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]); 
  free(ak2);
  free(ak3);
  free(ak4);
  free(ak5);
  free(ak6);
  free(ytemp);
  //  free_vector(ytemp,1,n); 
  //  free_vector(ak6,1,n); 
  //  free_vector(ak5,1,n); 
  //  free_vector(ak4,1,n); 
  //  free_vector(ak3,1,n); 
  //  free_vector(ak2,1,n); 
} 
#undef NRANSI
