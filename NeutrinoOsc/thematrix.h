/* For the time being matrices are defined as a 3 x 3 array of complex 
numbers*/
/* A matrix will be 

           j=1 j=2 j=3

i=1        a1  a2  a3 
i=2        b1  b2  b3 
i=3        c1  c2  c3 
*/

#include "complex.h"

#define PRECISION 1.e-8  /* Precision to which elements are added in a power 
			    expansion  */

#define ROWS    3
#define COLUMNS 3

typedef struct MATRIX matrix;
struct MATRIX {
  complex m[ROWS][COLUMNS];
};

void mcreate(matrix *A,complex c[COLUMNS][ROWS]);
void midentity(matrix *A );
void mreset(matrix *A);
matrix madd(matrix A, matrix B);
matrix mdiff(matrix A, matrix B);
matrix mscale(matrix A, complex c);
matrix mmul(matrix A, matrix B);
matrix mpow(matrix A, int n);
matrix mstar(matrix A);
matrix mtrans(matrix A);
matrix mdagger(matrix A);
matrix mexp(matrix A);
matrix mcos(matrix A);
matrix msin(matrix A);
matrix mmod_of_elements(matrix A);
double mmax_element(matrix A);
double mget_prob(matrix A, int i, int j);
void mprint(matrix A);
double fact(int n);
