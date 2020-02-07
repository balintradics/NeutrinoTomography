/**************************************************************************/
/**************************************************************************/
/* 
             UTILITY PROGRAM THAT ALLOWS TO PERFORM SOME SIMPLE 
                       ALGEBRA WITH COMPLEX MATRICES
 */
/**************************************************************************/
/**************************************************************************/
#include "stdio.h"
#include "math.h"

#include "thematrix.h"

/*\docF*/
void mcreate(matrix *A,complex c[ROWS][COLUMNS])
/* Creates for the time being matrices defined as a 3 x 3 array of complex 
numbers. A matrix will be 

           j=1 j=2 j=3

i=1        a1  a2  a3 
i=2        b1  b2  b3 
i=3        c1  c2  c3 
*/
{
  int i,j;

  for(i=0;i<ROWS;i++)
    {
      for(j=0;j<COLUMNS;j++)
	{
	  A->m[i][j].Re =c[i][j].Re;
	  A->m[i][j].Im =c[i][j].Im;
	}
    }
}

/*\docF*/
void midentity(matrix *A) 
/*Creates the identity matrix */
{
  int i,j;
  
  for(i=0;i<ROWS;i++)
    for(j=0;j<COLUMNS;j++)
      A->m[i][j].Re =A->m[i][j].Im = 0.;

  A->m[0][0].Re= A->m[1][1].Re= A->m[2][2].Re = 1.;
}

/*\docF*/
void mreset(matrix *A) 
/*Reset matrix elements */
{
  int i,j;
  
  for(i=0;i<ROWS;i++)
    for(j=0;j<COLUMNS;j++)
      A->m[i][j].Re =A->m[i][j].Im = 0.;
}

/*\docF*/
matrix madd(matrix A, matrix B)
/*Adds two matrices */
{
  int i,j;
  matrix C;

  for(i=0;i<ROWS;i++)
    for(j=0;j<COLUMNS;j++)
      C.m[i][j] = cadd(A.m[i][j], B.m[i][j]);
  return C;
}

/*\docF*/
matrix mdiff(matrix A, matrix B)
/*Subtract two matrices */
{
  int i,j;
  matrix C;

  for(i=0;i<ROWS;i++)
    for(j=0;j<COLUMNS;j++)
      C.m[i][j] = cdiff(A.m[i][j], B.m[i][j]);
  return C;
}

/*\docF*/
matrix mscale(matrix A, complex c)
/*Multiply a matrix by a complex number */
{
  int i,j;
  matrix B;
  for(i=0;i<ROWS;i++)
    for(j=0;j<COLUMNS;j++)
      B.m[i][j] = cmul(A.m[i][j],c);
  
  return B;
}

/*\docF*/
matrix mstar(matrix A)
/*Compute the matrix complex conjugate */
{
  matrix Astar;
  int i,j;

  for(i=0;i<ROWS;i++)
    for(j=0;j<COLUMNS;j++)
      Astar.m[i][j] = cconj(A.m[i][j]);

  return Astar;
}

/*\docF*/
matrix mtrans(matrix A)
/* Transpose the matrix */
{
  int i,j;
  matrix B;

  for(i=0;i<ROWS;i++)
    for(j=0;j<COLUMNS;j++)
      B.m[j][i] = A.m[i][j];

  return B;
}

/*\docF*/
matrix mdagger(matrix A)
/* Computes transpose complex conjugate of a matrix */
{
  matrix B;
  
  B = mstar(mtrans(A));
  return B;
}

/*\docF*/
matrix mmul(matrix A, matrix B) 
/* Multipliy two matrices */
{
  int i,j,k;
  matrix C;

  mreset(&C);

  for(i=0;i<ROWS;i++)      
    for(j=0;j<COLUMNS;j++)	  
      for(k=0;k<ROWS;k++)	      
	C.m[i][j]= cadd(C.m[i][j],(cmul(A.m[i][k], B.m[k][j])));

  return C;
}

/*\docF*/
matrix mpow(matrix A, int n) 
/*Computes A to the n */
{
  matrix B,C;

  B=A;

  if(!n) 
    printf("No way to do this operation\n");

  while(--n) 
    {
      C=mmul(B,A);
      B=C;
    }

  return B;
}

/*\docF*/
matrix mexp(matrix A) 
/* Computes the exponential of a matrix through a power expansion serie */
{
  /* We have to write the power expansion of the exponential 
     e^A = I + A + (1/2!) A^2 + (1/3!) A^3 + ...*/
  matrix B,C,D;
  int n;
  complex c;
  
  midentity(&C);

  /* First two terms of the serie */
  C=madd(C,A);
  B=A;

  /* Next terms of the expansion */
  for(n=2; ;n++)
    {
      D=mmul(B,A);
      B=D;
      
      c.Re=1./fact(n);
      c.Im=0.;
      D=mscale(D,c);

      if(mmax_element(D) <= PRECISION) 
	break;

      C=madd(C,D);
    }

  return C;
}
/*\docF*/
matrix mcos(matrix A) 
/* Computes the cosine of a matrix through a power expansion serie */
{
  /* We have to write the power expansion of the cosine
     cos(A) = I - (1/2!) A^2 + (1/4!) A^4 + ...*/
  matrix B,C,D,A2;
  int n,l;
  complex c;

  midentity(&C);

  l=1;
  B=C;
  A2=mmul(A,A);

  for(n=2; ;n+=2)
    {
      D=mmul(A2,B);
      B=D;
      
      c.Re=pow(-1,l)/fact(n);
      c.Im=0.;
      l++;

      D=mscale(D,c);
      
      if(mmax_element(D) <= PRECISION) 
	break;

      C=madd(C,D);
    }

  return C;
}

/*\docF*/
matrix msin(matrix A)
/* Computes the sine of a matrix through a power expansion serie */
{
  /* We have to write the power expansion of the cosine
     sin(A) = A - (1/3!) A^3 + (1/5!) A^5 + ...*/
  matrix A2,B,C,D;
  int n,l;
  complex c;
  
  /* First term */
  C=A;

  A2=mmul(A,A);
  B=A;

  /* Next terms */
  l=1;
  for(n=3; ;n+=2)
    {
      D=mmul(A2,B);
      B=D;
      
      c.Re=pow(-1,l)/fact(n);
      c.Im=0.;
      l++;

      D=mscale(D,c);

      if(mmax_element(D) <= PRECISION) 
	break;

      C=madd(C,D);
    }

  return C;
}

/*\docF*/
matrix mmod_of_elements(matrix A) 
/* Compute the modulus of each of the nine matrix elements */
{
  /* This returns a REAL matrix where the elements are the modulus 
     of the elements of the input imaginary matrix. 
     Useful to compute oscillation probabilities!! */
  matrix B;
  int i,j;

  mreset(&B);

  for(i=0;i<ROWS;i++)      
    for(j=0;j<COLUMNS;j++)	  
	B.m[i][j].Re=cmod(A.m[i][j]);

  return B;
}

double mget_prob(matrix A, int k, int l) 
{
  int i,j;
  double prob;

  /* The matrix should be real  otherwise we don't know what the hell to do!*/
  for(i=0;i<ROWS;i++)      
    for(j=0;j<COLUMNS;j++)	  
      if(abs(A.m[i][j].Im) != 0) 
	{
	  printf("mget_prob: Not a real matrix!!\n");
	  return -888.;
	}

  prob=A.m[k][l].Re *A.m[k][l].Re;
  return prob;
}

/*\docF*/
double mmax_element(matrix A) 
/* Returns the maximum of the moduli of the matrix elements */
{
  double modulus;
  double max=-999.;
  int i,j;

  for(i=0;i<ROWS;i++)      
    for(j=0;j<COLUMNS;j++)	  
      {
	modulus=cmod(A.m[i][j]);
	if(modulus > max) 
	  max = modulus;
      }
  return max;
}

/*\docF*/
void mprint(matrix A) 
/* Prints out the values of the matrix elements */
{
  int i,j;
  
  for(i=0;i<ROWS;i++)
    {
      printf("Print row %3d\n",i);
      for(j=0;j<COLUMNS;j++)
	{
	  printf("%.20e + i %.20e   ",A.m[i][j].Re,A.m[i][j].Im);
	}
      printf("  \n");
    }
  printf("  \n");
}

/* Certainly not matrix algebra, but useful to compute some nice stuff */
double fact(int n) 
{
  double f=1.;

  while(n) 
      f*=n--;

  return f;
}
