/**************************************************************************/
/**************************************************************************/
/* 
             UTILITY PROGRAM THAT ALLOWS TO PERFORM SOME SIMPLE 
                       ALGEBRA WITH COMPLEX NUMBERS
 */
/**************************************************************************/
/**************************************************************************/
#include "stdio.h"
#include "math.h"

#include "complex.h"

complex I;

/*\docF*/
void init_cmath(void)
/* Initializes real part to zero and imaginary part to one */
{
  I.Re=0;
  I.Im=1;
}

/*\docF*/
complex cadd(complex a, complex b)
/*Adds two complex numbers */
{
  complex c;
  CADDAB(c,a,b);
  return c;
}

/*\docF*/
complex cdiff(complex a, complex b)
/*Subtract two complex numbers */
{
  complex c;
  c.Re = a.Re-b.Re;
  c.Im = a.Im-b.Im;
  return c;
}

/*\docF*/
complex cmul(complex a, complex b)
/* Multiplication of two complex numbers */
{
  complex c;
  CMULAB(c,a,b);
  return c;
}

/*\docF*/
complex cmulre(double a, complex b)
/* Multiplication of real and complex numbers */
{
  complex c;
  c.Re=a*b.Re;
  c.Im=a*b.Im;
  return c;
}

/*\docF*/
complex cexp(complex a)
/*Exponential of a complex number */
{
  complex v;
  v.Re=exp(a.Re)*cos(a.Im);
  v.Im=exp(a.Re)*sin(a.Im);
  return v;
}

/*\docF*/
double cmod(complex a)
/*Modulus of a complex number */
{
  double m;
  m=CMOD(a);
  return m;
}

/*\docF*/
complex cconj(complex a) 
/*Complex conjugate */
{
  complex v;
  v.Re = a.Re; 
  v.Im = -a.Im;
  return v;
}

/*\docF*/
void creset(complex* a) 
/*Complex conjugate */
{
  a->Re = 0.;
  a->Im = 0.;
}
