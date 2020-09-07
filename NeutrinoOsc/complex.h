/*  */

#ifndef _COMPLEX_H_
#define _COMPLEX_H_

#include "math.h"

#define CMULAB(C,A,B)   C.Re = A.Re*B.Re - A.Im*B.Im;  C.Im = A.Re*B.Im + A.Im*B.Re;
#define CMULACJB(C,A,B) C.Re = A.Re*B.Re + A.Im*B.Im;  C.Im = - A.Re*B.Im + A.Im*B.Re;
#define CADDAB(C,A,B)   C.Re = A.Re+B.Re; C.Im = A.Im+B.Im;
#define CMULAB(C,A,B)   C.Re = A.Re*B.Re - A.Im*B.Im;  C.Im = A.Re*B.Im + A.Im*B.Re;
#define CMOD(A)         sqrt(A.Re*A.Re+A.Im*A.Im)
#define CMOD2(A)        (A.Re*A.Re+A.Im*A.Im)
#define CEXP(A)         {double tempcexp; tempcexp=exp(A.Re);A.Re=tempcexp*cos(A.Im);A.Im=tempcexp*sin(A.Im);}
#define CEXPIMAG(A)     A.Re=cos(A.Im);A.Im=sin(A.Im);


typedef struct COMPLEX complex;
struct COMPLEX {
  double Re;
  double Im;
};

extern complex I;
void init_cmath(void);
complex cadd(complex a, complex b);
complex cdiff(complex a, complex b);
complex cmul(complex a, complex b);
complex cmulre(double a, complex b);
complex cexp(complex a);
double cmod(complex a);
complex cconj(complex a); 
void creset(complex *a); 

#endif
