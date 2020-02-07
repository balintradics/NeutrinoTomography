/* Neutrino oscillations */
/* The "most important code" written in 1999. */
/* You give us the matrix, we tell you if it's excluded by data! */
/* by ARABs */
/* Send comments to A. Rubbia */

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define SQR(a) ((a)*(a))

#define CUBE(a) ((a)*(a)*(a))

#include "neutrino_osc.h"
#include "odeint.h"


// declarations moved to the cpp file
// to avoid duplicate symbols
int Recompute_Hamiltonian;
NeutrinoOsc NeuOsc;
NeutrinoOsc MatOsc;

// Just a comment 
double hbarc = 197.327; /* MeV * fm */

static int init_parms = 1;

static int awake_the_master;


//Following stuff needed for variable density integration
matrix Hmatter;


int kount = 0;
double *xp,**yp,dxsav; /* defined in odeint.c */ 

void rkqs(float y[], float dydx[], int n, float *x, float htry, float eps, float yscal[], float *hdid, float *hnext, void (*derivs)(float, float [], float []));

/*************************************************************************************/
/*                    P A R A M E T E R S       S E T T I N G                        */
/*************************************************************************************/

void nuox_dump_mix_matrix(char *mess, complex U[NDIM][NDIM])
{
  int alpha,i;
  printf("%s\n",mess);
  for(alpha=0;alpha<3;alpha++)
    {
      for(i=0;i<3;i++)
	printf("%10.3e  ",CMOD2(U[alpha][i]));
      printf("\n");
    }
}

/*\docF*/
void nuox_dump_a_osc_matrix(NeutrinoOsc *a)
/* Dump an oscillation matrix */
{
  printf("    Neutrino masses:   dm21=%10.3e   dm32=%10.3e\n",
	 a->dm21,a->dm32);
  printf("    Oscillation angles: t12=%10.3e  t13=%10.3e  t23=%10.3e     delta=%10.3e\n",
	 a->t12,a->t13,a->t23,a->delta);
  nuox_dump_mix_matrix((char*)("|U|**2:"),NeuOsc.U);
}


/*\docF*/
void nuox_dump_osc_matrix()
/* Dump the oscillation matrix */
{
  nuox_dump_a_osc_matrix(&NeuOsc);
}

/*\docF*/
void nuox_input_matrix_ckm_(double* dm32,double *dm21, double *t12,double *t13,double *t23,double *delta13)
{
nuox_input_matrix_CKM(*dm32,*dm21,*t12,*t13,*t23,*delta13);

}
void nuox_input_matrix_CKM(double dm32,double dm21, double t12,double t13,double t23,double delta13)
/* Specifies the problem of the neutrino masses and mixing. The mixing matrix is
specified by means of the CKM parametrization.*/
{
  int i,j;
  double c13,s13;
  double c12,s12;
  double c23,s23;
  int a,k,l;

  if(init_parms)
    {
      init_cmath();
      init_parms = 0;
    }

  NeuOsc.dm32=dm32;
  NeuOsc.dm21=dm21;
  NeuOsc.mass[0]=0;
  NeuOsc.mass[1]=sqrt(NeuOsc.mass[0]*NeuOsc.mass[0]+NeuOsc.dm21);
  NeuOsc.mass[2]=sqrt(NeuOsc.mass[1]*NeuOsc.mass[1]+NeuOsc.dm32);
  NeuOsc.dm31=NeuOsc.mass[2]*NeuOsc.mass[2]-NeuOsc.mass[0]*NeuOsc.mass[0];

  NeuOsc.t12 = t12;
  NeuOsc.t13 = t13;
  NeuOsc.t23 = t23;
  NeuOsc.delta = delta13;

  NeuOsc.sin22t12 = SQR(sin(2.*t12));
  NeuOsc.sin22t13 = SQR(sin(2.*t13));
  NeuOsc.sin22t23 = SQR(sin(2.*t23));

  NeuOsc.sin2t12 = SQR(sin(t12));
  NeuOsc.sin2t13 = SQR(sin(t13));
  NeuOsc.sin2t23 = SQR(sin(t23));

  c13 = cos(NeuOsc.t13);
  s13 = sin(NeuOsc.t13);
  
  c12 = cos(NeuOsc.t12);
  s12 = sin(NeuOsc.t12);
  
  c23 = cos(NeuOsc.t23);
  s23 = sin(NeuOsc.t23);
  
  for(i=0;i<NDIM;i++)
    for(j=0;j<NDIM;j++)
      NeuOsc.U[i][j].Re=NeuOsc.U[i][j].Im=0;

  //All the elements of the mixing matrix are defined here
  // Real part
  NeuOsc.U[0][0].Re= c12*c13;                                     /* e1 */
  NeuOsc.U[0][1].Re= s12*c13;                                     /* e2 */
  NeuOsc.U[0][2].Re= s13*cos(delta13);                            /* e3 */

  NeuOsc.U[1][0].Re= -s12*c23-(c12*s23*s13*cos(delta13));         /* mNeuOsc.U1 */
  NeuOsc.U[1][1].Re= c12*c23-(s12*s23*s13*cos(delta13));          /* mNeuOsc.U2 */
  NeuOsc.U[1][2].Re= s23*c13;                                     /* mNeuOsc.U3 */

  NeuOsc.U[2][0].Re= s12*s23-(c12*c23*s13*cos(delta13));          /* taNeuOsc.U1 */
  NeuOsc.U[2][1].Re= -c12*s23-(s12*c23*s13*cos(delta13));         /* taNeuOsc.U2 */
  NeuOsc.U[2][2].Re= c23*c13;                                     /* taNeuOsc.U3 */

  // Imaginary part
  NeuOsc.U[0][0].Im= 0.;                                           /* e1 */
  NeuOsc.U[0][1].Im= 0.;                                           /* e2 */
  NeuOsc.U[0][2].Im= -s13*sin(delta13);                            /* e3 */

  NeuOsc.U[1][0].Im= -c12*s23*s13*sin(delta13);                    /* mNeuOsc.U1 */
  NeuOsc.U[1][1].Im= -s12*s23*s13*sin(delta13);                    /* mNeuOsc.U2 */
  NeuOsc.U[1][2].Im= 0.;                                           /* mNeuOsc.U3 */

  NeuOsc.U[2][0].Im= -c12*c23*s13*sin(delta13);                   /* taNeuOsc.U1 */
  NeuOsc.U[2][1].Im= -s12*c23*s13*sin(delta13);                   /* taNeuOsc.U2 */
  NeuOsc.U[2][2].Im= 0.;                                           /* tau3 */


  NeuOsc.mass2[0]=NeuOsc.mass[0]*NeuOsc.mass[0];
  NeuOsc.mass2[1]=NeuOsc.mass[1]*NeuOsc.mass[1];
  NeuOsc.mass2[2]=NeuOsc.mass[2]*NeuOsc.mass[2];
  for(a=0;a<NDIM;a++)
    for(k=0;k<NDIM;k++)
      for(l=0;l<NDIM;l++)
	{      
	  CMULACJB(NeuOsc.Ukla[a][k][l],NeuOsc.U[k][a],NeuOsc.U[l][a]);
	}
}

/*\docF*/
void nux_CKM_from_U(NeutrinoOsc *a)
/* Fill CKM parametrization from U elements */
{
  double s132,s122,s232;
  s132 = CMOD2(a->U[0][2]);
  s122 = CMOD2(a->U[0][1])/(1.0-s132);
  s232 = CMOD2(a->U[1][2])/(1.0-s132);
  a->t13=asin(sqrt(s132));
  a->t12=asin(sqrt(s122));
  a->t23=asin(sqrt(s232));

  a->sin22t12 = SQR(sin(2.*a->t12));
  a->sin22t13 = SQR(sin(2.*a->t13));
  a->sin22t23 = SQR(sin(2.*a->t23));

  a->sin2t12 = SQR(sin(a->t12));
  a->sin2t13 = SQR(sin(a->t13));
  a->sin2t23 = SQR(sin(a->t23));
}

/*\docF*/
void nuox_set_new_CKM_matrix(double m1_2,double m2_2, double m3_2, double t12,double t13,double t23,double delta13)
/* Redefine mass differences and mixing angles once mass eigenvalues and 
new mixing angles have been found for propagation in matter of constant density */
{
  int i,j;
  double c13,s13;
  double c12,s12;
  double c23,s23;
  int a,k,l;

  MatOsc = NeuOsc;

  /* And now replace with the new computed values for mass eigenstates and 
     mixing angles */
  MatOsc.mass2[0]=m1_2;
  MatOsc.mass2[1]=m2_2;
  MatOsc.mass2[2]=m3_2;

  MatOsc.mass[0]=sqrt(MatOsc.mass2[0]);
  MatOsc.mass[1]=sqrt(MatOsc.mass2[1]);
  MatOsc.mass[2]=sqrt(MatOsc.mass2[2]);

  MatOsc.dm31= MatOsc.mass2[2]- MatOsc.mass2[0];
  MatOsc.dm21= MatOsc.mass2[1]- MatOsc.mass2[0];
  MatOsc.dm32= MatOsc.mass2[2]- MatOsc.mass2[1];

  MatOsc.t12 = t12;
  MatOsc.t13 = t13;
  MatOsc.t23 = t23;
  MatOsc.delta = delta13;

  MatOsc.sin22t12 = SQR(sin(2.*t12));
  MatOsc.sin22t13 = SQR(sin(2.*t13));
  MatOsc.sin22t23 = SQR(sin(2.*t23));

  MatOsc.sin2t12 = SQR(sin(t12));
  MatOsc.sin2t13 = SQR(sin(t13));
  MatOsc.sin2t23 = SQR(sin(t23));

  c13 = cos(MatOsc.t13);
  s13 = sin(MatOsc.t13);
  
  c12 = cos(MatOsc.t12);
  s12 = sin(MatOsc.t12);
  
  c23 = cos(MatOsc.t23);
  s23 = sin(MatOsc.t23);
  
  for(i=0;i<NDIM;i++)
    for(j=0;j<NDIM;j++)
      MatOsc.U[i][j].Re=MatOsc.U[i][j].Im=0;

  //All the elements of the mixing matrix are defined here
  // Real part
  MatOsc.U[0][0].Re= c12*c13;                                     /* e1 */
  MatOsc.U[0][1].Re= s12*c13;                                     /* e2 */
  MatOsc.U[0][2].Re= s13*cos(delta13);                            /* e3 */

  MatOsc.U[1][0].Re= -s12*c23-(c12*s23*s13*cos(delta13));         /* mMatOsc.U1 */
  MatOsc.U[1][1].Re= c12*c23-(s12*s23*s13*cos(delta13));          /* mMatOsc.U2 */
  MatOsc.U[1][2].Re= s23*c13;                                     /* mMatOsc.U3 */

  MatOsc.U[2][0].Re= s12*s23-(c12*c23*s13*cos(delta13));          /* taMatOsc.U1 */
  MatOsc.U[2][1].Re= -c12*s23-(s12*c23*s13*cos(delta13));         /* taMatOsc.U2 */
  MatOsc.U[2][2].Re= c23*c13;                                     /* taMatOsc.U3 */

  // Imaginary part
  MatOsc.U[0][0].Im= 0.;                                           /* e1 */
  MatOsc.U[0][1].Im= 0.;                                           /* e2 */
  MatOsc.U[0][2].Im= -s13*sin(delta13);                            /* e3 */

  MatOsc.U[1][0].Im= -c12*s23*s13*sin(delta13);                    /* mMatOsc.U1 */
  MatOsc.U[1][1].Im= -s12*s23*s13*sin(delta13);                    /* mMatOsc.U2 */
  MatOsc.U[1][2].Im= 0.;                                           /* mMatOsc.U3 */

  MatOsc.U[2][0].Im= -c12*c23*s13*sin(delta13);                   /* taMatOsc.U1 */
  MatOsc.U[2][1].Im= -s12*c23*s13*sin(delta13);                   /* taMatOsc.U2 */
  MatOsc.U[2][2].Im= 0.;                                           /* tau3 */

  for(a=0;a<NDIM;a++)
    for(k=0;k<NDIM;k++)
      for(l=0;l<NDIM;l++)
	{      
	  CMULACJB(MatOsc.Ukla[a][k][l],MatOsc.U[k][a],MatOsc.U[l][a]);
	}
}

/*\docF*/
void nuox_input_matrix_CKM_matter(int i)
/* Specifies the problem of the neutrino masses and mixing. The mixing matrix is
specified by means of the CKM parametrization. Propagation in matter with constant density assumed */
{
  double A,B,C,D,S,angle;
  double A2,S2,E2,F2;
  double factor,factor1,factor2,factor3;
  double alpha,beta,E,F;
  double DM32,DM21,DM31,t12m,t13m,t23m,delta13m;
  double M1_2,M2_2,M3_2;

  double c13,s13,c132,s132;
  double c12,s12,c122,s122;
  double c23,s23,c232,s232;
  
  double sm12_2,sm13_2,sm23_2;
  double cdeltam,sdeltam;
  double BB;
  double denom;

  int iflag;

  /* I have to compute the hamiltonian eigenvalues and mixing angles in 
     matter 
     Formulae are extracted from H.W. Zaglauer and K.H. Schwarzer 
                                 Z. Phys C 40 (1988) 273   */

  /* Useful values from vacuum mixing angles */
  c13 = cos(NeuOsc.t13);
  s13 = sin(NeuOsc.t13);
  c132 = SQR(c13);
  s132 = SQR(s13);

  c12 = cos(NeuOsc.t12);
  s12 = sin(NeuOsc.t12);
  c122 = SQR(c12);
  s122 = SQR(s12);

  c23 = cos(NeuOsc.t23);
  s23 = sin(NeuOsc.t23);
  c232 = SQR(c23);
  s232 = SQR(s23);

  /*********************************************************************/
  /*                   Effectives masses in matter                     */
  /*********************************************************************/

  /* This is the term: 2* sqrt(2) * G_F * n_e * E  in eV^2 */
  D = 7.56e-5*NeuOsc.density*NeuOsc.Energy;
  if(NeuOsc.neutype<0) /* Anti-neutrinos */
    D = -D;

  A = NeuOsc.dm21 + NeuOsc.dm31 + D;
  A2 = SQR(A);
  
  B = NeuOsc.dm21*NeuOsc.dm31+
    D*(NeuOsc.dm31*c132+NeuOsc.dm21*(c132*c122+s132));

  C = D*NeuOsc.dm21*NeuOsc.dm31*c132*c122;

  angle = (2.*CUBE(A)-9.*A*B+27.*C)/(2.*pow(A2-3.0*B,3./2.));
#ifdef DEBUG
  printf("angle %f\n",angle);
#endif
  angle=acos(angle);
  S = cos((1./3.)*angle);
  S2=SQR(S);

  /* The mass eigenvalues squared are...*/
  factor = NeuOsc.mass2[0] + (A/3.);
  factor1 = sqrt(A2-3*B)*S/3.;
  factor2 = sqrt(A2-3*B)*sqrt(1-S2)/3.;
  M1_2 = factor - factor1-sqrt(3.)*factor2;
  M2_2 = factor - factor1+sqrt(3.)*factor2;
  M3_2 = factor + 2*factor1;

  //#ifdef TEST
  //  if(M1_2<0)
  //    M1_2=0.;
  //#endif

#ifdef DEBUG
  printf("DD D=%10.4f   M1_2: %10.4e  M2_2: %10.4e  M3_2: %10.4e\n",
	 D,M1_2,M2_2,M3_2);
#endif

  /* Mass differences in matter */
  DM21=M2_2-M1_2;
  DM31=M3_2-M1_2;
  DM32=M3_2-M2_2;
#ifdef DEBUG
  printf("DD E=%10.4f   M21: %10.4e  M31: %10.4e  M32: %10.4e\n",
	 NeuOsc.Energy,DM21,DM31,DM32);
#endif

  /*********************************************************************/
  /* Now we have to define the angles in matter                        */
  /*********************************************************************/

  /* More useful definitions */
  alpha = NeuOsc.mass2[2]*c132+
          NeuOsc.mass2[1]*(c132*c122+s132)+
          NeuOsc.mass2[0]*(c132*s122+s132);

  beta = NeuOsc.mass2[2]*c132*
                  (NeuOsc.mass2[1]*c122+NeuOsc.mass2[0]*s122)+
         NeuOsc.mass2[1]*NeuOsc.mass2[0]*s132;

  E=c13*s13*(NeuOsc.dm31*(M3_2-NeuOsc.mass2[0]-NeuOsc.dm21) - 
             NeuOsc.dm21*(M3_2-NeuOsc.mass2[0]-NeuOsc.dm31)*s122);
  E2 = SQR(E);

  F = (M3_2-NeuOsc.mass2[0]-NeuOsc.dm31)*NeuOsc.dm21*c12*s12*c13;
  F2 = SQR(F);

  /* Now the squared sine of the mixing angles */
  sm12_2 = -(SQR(M2_2)-alpha*M2_2+beta)*DM31;
  sm12_2 = sm12_2/(DM32*(SQR(M1_2)-alpha*M1_2+beta)-
		   DM31*(SQR(M2_2)-alpha*M2_2+beta));

  sm13_2 = (SQR(M3_2)-alpha*M3_2+beta)/DM31/DM32;

  sm23_2 = E2*s232+F2*c232+2*E*F*c23*s23*cos(NeuOsc.delta);
  sm23_2 = sm23_2/(E2+F2);
  
  /* To avoid stupid crashes due to precision */
  iflag = 0;
  if(sm12_2 > 1.) 
    {
      sm12_2 = 1.;
      iflag =1;
    }
  if(sm13_2 > 1.) 
    {
      sm13_2 = 1.;
      iflag =1;
    }
  if(sm23_2 > 1.) 
    {
      sm23_2 = 1.;
      iflag =1;
    }
  
  if(sm12_2 < 0) 
    {
      sm12_2 = 0;
      iflag =1;
    }
  if(sm13_2 < 0) 
    {
      sm13_2 = 0;
      iflag =1;
    }
  if(sm23_2 < 0) 
    {
      sm23_2 = 0;
      iflag =1;
    }

#ifdef DEBUG
  if(iflag)
    printf("AA sm12_2 %10.4e sm13_2 %10.4e sm23_2 %10.4e \n",
	   sm12_2,sm13_2,sm23_2);
#endif

  t12m = asin(sqrt(sm12_2));
  t13m = asin(sqrt(sm13_2));
  t23m = asin(sqrt(sm23_2));

#ifdef DEBUG
  printf("AA t12m %10.4e(%10.4e)  t13m %10.4e(%10.4e)  t23m %10.4e(%10.4e)\n",
	 t12m,NeuOsc.t12,t13m,NeuOsc.t13,t23m,NeuOsc.t23);
#endif
  
  /*********************************************************************/
  /* And now the CP violation phase in matter !!!!*/
  /*********************************************************************/
  factor3 = 2*E*F*c23*s23*cos(NeuOsc.delta);
  cdeltam = (E2-F2)*cos(NeuOsc.delta)*s23*c23+E*F*(c232-s232);
  sdeltam = -(E2+F2)*sin(NeuOsc.delta)*s23*c23;
  denom = (E2*s232+F2*c232+factor3)*(E2*c232+F2*s232-factor3);
  if(denom<0)
    {
      printf("denom is negative??? %f\n",denom);
    }
  denom=sqrt(denom);

  if(denom>0)
    {
      cdeltam = cdeltam/denom;
#ifdef DEBUG
      if(fabs(cdeltam)>1.)
	{
	  printf("cdeltam to big! %f\n",cdeltam);
	}
#endif
      if(cdeltam>1.0)cdeltam=1.0;
      if(cdeltam<-1.0)cdeltam=-1.0;
      if(sdeltam*NeuOsc.neutype<0)
      	delta13m = acos(cdeltam);
      else
	delta13m = - acos(cdeltam);
    }
  else
    delta13m=NeuOsc.delta;


  if(NeuOsc.neutype<0) /* Anti-neutrinos */
    delta13m = - delta13m;

  /* Now we redifine the parameters governing the oscillation in 
     matter with constant density */
  nuox_set_new_CKM_matrix(M1_2,M2_2,M3_2,t12m,t13m,t23m,delta13m);
  //  nuox_set_new_CKM_matrix(NeuOsc.mass2[0],NeuOsc.mass2[1],NeuOsc.mass2[2],t12m,t13m,t23m,delta13m);
  //  nuox_dump_a_osc_matrix(&MatOsc);
}
  
void nuox_input_matrix_CKM_xing(void)
/* Specifies the problem of the neutrino masses and mixing. The mixing matrix is
specified by means of the CKM parametrization. Propagation in matter with constant density assumed */
{
  complex Um[NDIM][NDIM];

  double Ni[NDIM];
  double Di[NDIM];
  int i,j,k;
  int idxj[3],idxk[3];

  double lambda[3];

  int alpha;
  complex v1,v2,v3,v4;

  double A,A2,B,C,D,angle,S,S2;
  double factor,factor1,factor2;

  MatOsc = NeuOsc;

  D = 7.56e-5*NeuOsc.density*NeuOsc.Energy;
  if(NeuOsc.neutype<0) /* Anti-neutrinos */
    D = -D;

  A = NeuOsc.dm21 + NeuOsc.dm31 + D;
  A2 = SQR(A);
  
  B = NeuOsc.dm21*NeuOsc.dm31+
    D*(NeuOsc.dm31*(1-CMOD2(NeuOsc.U[0][2]))+NeuOsc.dm21*(1-CMOD2(NeuOsc.U[0][1])));

  C = D*NeuOsc.dm21*NeuOsc.dm31*CMOD2(NeuOsc.U[0][0]);

  angle = (2.*CUBE(A)-9.*A*B+27.*C)/(2.*pow(A2-3.0*B,3./2.));
  angle=acos(angle);
  S = cos((1./3.)*angle);
  S2=SQR(S);

  /* The mass eigenvalues squared are...*/
  factor = NeuOsc.mass2[0] + (A/3.);
  factor1 = sqrt(A2-3*B)*S/3.;
  factor2 = sqrt(A2-3*B)*sqrt(1-S2)/3.;
  MatOsc.mass2[0] = lambda[0] = factor - factor1-sqrt(3.)*factor2;
  MatOsc.mass2[1] = lambda[1] = factor - factor1+sqrt(3.)*factor2;
  MatOsc.mass2[2] = lambda[2] = factor + 2*factor1;

  MatOsc.mass[0]=sqrt(MatOsc.mass2[0]);
  MatOsc.mass[1]=sqrt(MatOsc.mass2[1]);
  MatOsc.mass[2]=sqrt(MatOsc.mass2[2]);

  MatOsc.dm31= MatOsc.mass2[2]- MatOsc.mass2[0];
  MatOsc.dm21= MatOsc.mass2[1]- MatOsc.mass2[0];
  MatOsc.dm32= MatOsc.mass2[2]- MatOsc.mass2[1];

  idxj[0]=1;idxj[1]=2;idxj[2]=0;
  idxk[0]=2;idxk[1]=0;idxk[2]=1;
  for(i=0;i<3;i++)
    {
      j=idxj[i];k=idxk[i];
      Ni[i]=(lambda[i]-NeuOsc.mass2[j])*(lambda[i]-NeuOsc.mass2[k])-D*((lambda[i]-NeuOsc.mass2[j])*CMOD2(NeuOsc.U[0][k])
							+(lambda[i]-NeuOsc.mass2[k])*CMOD2(NeuOsc.U[0][j]));
      Di[i]=SQR(Ni[i])+SQR(D)*CMOD2(NeuOsc.U[0][i])*(SQR(lambda[i]-NeuOsc.mass2[j])*CMOD2(NeuOsc.U[0][k])+
					    SQR(lambda[i]-NeuOsc.mass2[k])*CMOD2(NeuOsc.U[0][j]));
      Di[i]=sqrt(Di[i]);
    }

  for(alpha=0;alpha<3;alpha++)
    for(i=0;i<3;i++)
      {
	j=idxj[i];k=idxk[i];
	v1=cmulre(Ni[i]/Di[i],NeuOsc.U[alpha][i]);
	v2=cmul(cconj(NeuOsc.U[0][k]),NeuOsc.U[alpha][k]);
	v2=cmulre(lambda[i]-NeuOsc.mass2[j],v2);
	v3=cmul(cconj(NeuOsc.U[0][j]),NeuOsc.U[alpha][j]);
	v3=cmulre(lambda[i]-NeuOsc.mass2[k],v3);
	v4=cadd(v2,v3);
	v4=cmul(v4,NeuOsc.U[0][i]);
	v4=cmulre(D/Di[i],v4);
	MatOsc.U[alpha][i]=cadd(v1,v4);
      }
  nux_CKM_from_U(&MatOsc);
}

/*\docF*/
void nuox_set_propag_level_(int *ilevel, double *rho)
{
  nuox_set_propag_level(*ilevel, *rho);	
}

void nuox_set_propag_level(int ilevel, double rho)
/* Set the kind of propagation approximation 
 ilevel=0 : vacuum 
       =1 : matter, constant 
       =2 : matter, Earth profile approximated 
       =3 : matter, Earth profile exact (VERY SLOW) */
{
  /* Medium density */
  NeuOsc.density= rho;
  NeuOsc.ilevel = ilevel;
  if(NeuOsc.ilevel)
    nuox_set_hamiltonian(NeuOsc.ilevel);// was "1"
}

void nuox_set_hamiltonian(int i)
{
  /* If Recompute_Hamiltonian is equal to zero we compute oscillation 
     probabilities in matter without recomputing the hamiltonian 
     save some CPU */
  Recompute_Hamiltonian = i;
}

/*\docF*/
void nuox_set_neutrino_(double *l,double *e,int *ntyp)
{
nuox_set_neutrino(*l,*e,*ntyp);
}
void nuox_set_neutrino(double l,double e,int ntyp)

/* Specifies the path-length, the energy, the matter density and
   the type of neutrino. ntyp=1 for neutrinos, -1 for antineutrinos */
{
  // Length and energy 
  NeuOsc.Length= l;
  NeuOsc.Energy= e;  
  // Neutrino(=1) Anti-Neutrino(=-1)
  NeuOsc.neutype= ntyp;  

  /* Make sure for matter oscillations you recompute the hamiltonian 
     after changing neutrino parameters */
  if(NeuOsc.ilevel)
    nuox_set_hamiltonian(1);

  /* To save CPU time we have this flag */
  awake_the_master = 1;
}
void nuox_set_neutrino_costheta(double costheta) 
{
  NeuOsc.CosTheta=costheta;
}
/*\docF*/
static double nuox_propag(NeutrinoOsc *n, int l, int k)
/* This routine returns the probability for neutrino of type i 
   to be converted into neutrino type j after propagation in ``vacuum'' 
   0=electron, 1=muon,2=tau */  
{
  int a;
  complex expon;
  complex amplitude;
  double probability;
  complex c1;
  double factor;
  complex matele;

  amplitude.Re=0.;
  amplitude.Im=0.;
  factor = n->Length/(2e9*n->Energy)/(hbarc*1e6*1e-18);
  for(a=0;a<NDIM;a++)
    {      
      expon.Im=-factor*n->mass2[a];
      CEXPIMAG(expon);
      matele=n->Ukla[a][k][l];
      if(n->neutype<0)
	matele.Im=-matele.Im;
      //      printf("a %d   --> matele %10.4e %10.4e mod2 %10.4e\n",
      //	     a,matele.Re,matele.Im,CMOD2(matele));
      CMULAB(c1,matele,expon);
      CADDAB(amplitude,amplitude,c1);
      //      printf("a %d   --> amplitude %10.4e %10.4e mod2 %10.4e\n",
      //	     a,amplitude.Re,amplitude.Im,CMOD2(amplitude));
    }
  //  printf("amplitude %10.4e %10.4e \n",amplitude.Re,amplitude.Im);
  probability = CMOD2(amplitude);
  return probability;
}

/*\docF*/
matrix nuox_hamiltonian_matter(matrix Mmix,double mass[NDIM],double E,double L,int ipar) 
/*  In this subroutine we compute: -i*H*L/(hbar *c)
    where H is defined as ($U*H_masses*U^-1+H_matter)/2*E$
    ipar=1 means neutrino, ipar=-1 means anti-neutrino */
{
  complex c,p;
  matrix H,H1;
  matrix Mcos,Msin;
  matrix Mmass,MexpH,Minducedmass,Mdag;
  double factor_in_exp=-L/(hbarc*1e6*1e-18);
  
  // Create mass matrix
  mreset(&Mmass);
  Mmass.m[0][0].Re=mass[0]*mass[0];
  Mmass.m[1][1].Re=mass[1]*mass[1];
  Mmass.m[2][2].Re=mass[2]*mass[2];
  
  Mdag = mdagger(Mmix);
  H1=mmul(Mmass,Mdag);
  H=mmul(Mmix,H1);
  
  p.Im=0.;
  p.Re=1./2./NeuOsc.Energy/1.e9;
  H=mscale(H,p);
  
  c.Re=0.;
  c.Im=factor_in_exp; 
  H=mscale(H,c); 


  // Compute factor = 2 *sqrt(2) * G_F * Y_e * m_n^-1 * density * E
  // G_F = 8.92 e-8 GeV fm^3 
  // Y_e= 0.5 (this is an assumption!) 
  // m_n = 1.672 e-24 g 
  // density is given in g/cm^3 
  // energy is given in GeV
  //Create matrix associated to matter effects
  if(ipar<0) /* anti-neutrinos*/
    E=-E;
  mreset(&Minducedmass);
  Minducedmass.m[0][0].Re=0.;
  Minducedmass.m[0][0].Im= 7.577e-5*NeuOsc.density*E*factor_in_exp;
  Minducedmass=mscale(Minducedmass,p);

  //Sum the two contributions to the total hamiltonian 
  H=madd(H,Minducedmass);

  MexpH=mexp(H);
  
  return MexpH;
}

/*\docF*/
matrix nuox_amplitude_matter(void )
/* Neutrino propagation in matter using matrix formalism. It returns a matrix containing 
     the amplitudes for the nine possible oscillation processes */
{

  /**********************************************************************************/
  /**********************************************************************************/
  /* 
                 N E U T R I N O   P R O P A G A T I O N   I N   M A T T E R 
                             M A T R I X    F O R M A L I S M 
  */
  /**********************************************************************************/
  /**********************************************************************************/
  matrix Mmix,Hmat,Hmat_contri;
  int i;
  double DIV=1.;

  /* Create mixing matrix */
  mcreate(&Mmix,NeuOsc.U);
  
  /* Compute the Hamiltonian in matter */
  midentity(&Hmat);
  /* since we have problems when dealing with large distances, we split into smaller pieces*/
  for(i=0;i<DIV;i++)
    {
      Hmat_contri=nuox_hamiltonian_matter(Mmix,NeuOsc.mass,NeuOsc.Energy,NeuOsc.Length/DIV,NeuOsc.neutype);
      Hmat=mmul(Hmat,Hmat_contri);
    }

  return Hmat;
}


// ===============================================================================================
// ===============================================================================================
//                CODE TO COMPUTE OSCILLATIONS IN A MEDIUM OF VARYING DENSITTY
// ===============================================================================================
// ===============================================================================================
void locate(double xx[], int n, double x, int *j) 
     //Given an array xx[1..n], and given a value x, returns a value j such that x is between xx[j] and xx[j+1]. xx must be monotonic, either increasing or decreasing. j=0 or j=n is returned to indicate that x is out of range. 
{ 
  int ju,jm,jl; 
  int ascnd;

  jl=0; 
  //  Initialize lower 
  ju=n; 
  //and upper limits. 
  ascnd=(xx[n-1] >= xx[0]); 
  while (ju-jl > 1) 
    { 
      //  If we are not yet done, 
      jm=(ju+jl) >> 1; 
      // compute a midpoint, 
      if (x >= xx[jm] == ascnd) 
	jl=jm; //and replace either the lower limit 
      else 
	ju=jm; // or the upper limit, as appropriate. 
    } 
  //  Repeat until the test condition is satis ed. 
  if (x == xx[0]) 
    *j=0; 
  else if(x == xx[n-1]) //   Then set the output 
    *j=n-1; 
  else 
    *j=jl; 
}

double density (double x) 
{
  int j=0;

  double a,b;

  double depth[]={0.   ,15.  ,60.  ,100. ,200. ,300. ,350. ,400. ,413. ,500. ,
                  600. ,650. ,800. ,984. ,1000.,1200.,1400.,1600.,1800.,2000.,
                  2200.,2400.,2600.,2800.,2878.,3000.,3200.,3400.,3600.,3800.,
                  4000.,4200.,4400.,4600.,4800.,4982.,5000.,5121.,5200.,5400.,
                  5600.,5800.,6000.,6200.,6371.}; // in Km.

  double rho[]={2.840 ,2.840 ,3.332 ,3.348 ,3.387 ,3.424 ,3.441 ,3.775 ,3.795 ,3.925 ,
		4.075 ,4.150 ,4.380 ,4.529 ,4.538 ,4.655 ,4.768 ,4.877 ,4.983 ,5.087 ,
                5.188 ,5.288 ,5.387 ,5.487 ,5.527 ,10.121,10.421,10.697,10.948,11.176,
		11.383,11.570,11.737,11.887,12.017,12.121,12.130,12.197,12.760,12.855,
                12.916,13.02,13.059,13.08,13.09}; // g/cm^3. 

  int size=sizeof(depth)/sizeof(double);

  double Dist,R,R_E=6371.; // in Km
  double R_OuterCore = 3480.; // in Km



  //  ESTOY JODIENDO ALGO EN EL PUTO RETURN DE ESTA RUTINA!!!!!!!!!
  //  return 3.3;


  // Compute the distance to Earth's center to get the density
  //Dist=sqrt( SQR(NeuOsc.Length-x) + SQR(R_E) - 2.*R_E*(NeuOsc.Length-x)*fabs(NeuOsc.CosTheta) );
  // |(\vec{R_E} - (\vec{L - x}) )|

  // Simple geometrical distance from the Earth's surface
  Dist=R_E - NeuOsc.Length + x;
  if((Dist-R_E)>0.)
    return 0.; // Most likely propragation through the atmosphere: density value=0.

  // SLOW PROCEDURE: Read density values in a table and interpolate

  // If R=R_E this means we're in the crust and therefore depth=0, so...
  // the argument to locate is R_E-R and not simply R
  //  locate(depth,size,(R_E-R),&j);

  //  a=(rho[j+1]-rho[j])/(depth[j+1]-depth[j]);
  //  b=rho[j+1]-a*depth[j+1];

  //  return (a*x+b);


  // Parameterization of the Earth's density from I.Mocioui and R. Shrock
  //                                              Phys.Rev.D62:053017,2000
  R=Dist/R_E;

  double rho_ret = 2.6;

/**
  if(Dist>0. && Dist<1221.5) 
    rho_ret = (13.0885-8.8381*R*R);
  else if (Dist >=1221.5 && Dist<3480.)
    rho_ret = (12.5815-1.2638*R-3.6426*R*R-5.5281*R*R*R);
  else if (Dist>=3480. && Dist<5701.)
    rho_ret = (7.9565-6.4761*R + 5.5283*R*R-3.0807*R*R*R);
  else if (Dist >=5701.0 && Dist<5771.0)
    rho_ret = (5.3197-1.4836*R);
  else if (Dist>=5771.0 && Dist<5971.0)
    rho_ret = (11.2494-8.0298*R);
  else if (Dist>=5971.0 && Dist<6151.0)
    rho_ret = (7.1089-3.8045*R);
  else if (Dist>=6151.0 && Dist<6346.6)
    rho_ret = (2.6910+0.6924*R);
  else if (Dist>=6346.6 && Dist<6356.6)
    rho_ret = 2.9;
  else if (Dist>=6356.6 && Dist <=R_E)
    rho_ret = 2.6;
**/

  if(Dist>0. && Dist<1221.5) 
    rho_ret = (13.0885-8.8381*R*R);
  else if (Dist >=1221.5 && Dist<3480.)
    rho_ret = (12.5815-1.2638*R-3.6426*R*R-5.5281*R*R*R);
  else if (Dist>=3480. && Dist<5701.)
    rho_ret = (7.9565-6.4761*R + 5.5283*R*R-3.0807*R*R*R);
  else if (Dist >=5701.0 && Dist<5771.0)
    rho_ret = (5.3197-1.4836*R);
  else if (Dist>=5771.0 && Dist<5971.0)
    rho_ret = (11.2494-8.0298*R);
  else if (Dist>=5971.0 && Dist<6151.0)
    rho_ret = 19.3;
  else if (Dist>=6151.0 && Dist<6346.6)
    rho_ret = (2.6910+0.6924*R);
  else if (Dist>=6346.6 && Dist<6356.6)
    rho_ret = 2.9;
  else if (Dist>=6356.6 && Dist <=R_E)
    rho_ret = 2.6;

//  printf("CostTheta: %10.4e, NeutL: %10.4e, x: %10.4e, Dist: %10.4e, Dist-R_E: %10.4e, density: %10.4e\n",  NeuOsc.CosTheta, NeuOsc.Length, x, Dist, Dist-R_E, rho_ret);

  //default return 2.6
  return rho_ret;
}

void derivs(double x,double y[],double dydx[])
{
  int i,j,kk=0;
  double matterparam;
  // Compute factor = 2.*sqrt(2) * G_F * Y_e * rho * E /m_p
  // energy is given in GeV, distance in Km, mass in grams
  // to get right dimensions we compute parammsw = factor/(2 * hbarc * E)
  // and the value entering the calculation is...
  double parammsw=1.9207e-4;
  //  double parammsw=9.6035e-5;
 
  complex a[NDIM],b[NDIM];

  // Additional potential in matter
  matterparam= parammsw*density(x)*NeuOsc.neutype;
  Hmatter.m[0][0].Re+=matterparam;

  // Create complex numbers to be used in function definition
  for(j=0;j<NDIM;j++)
    {
      dydx[j]=dydx[j+NDIM]=0.;
      a[j].Re=y[j];
      a[j].Im=y[j+NDIM];
      creset(&b[j]);
    }

  // The function to be integrated
  for(i=0;i<NDIM;i++)
    {
      for(j=0;j<NDIM;j++)
	{
	  b[i]=cadd(b[i],cmul(Hmatter.m[i][j],a[j]));
	}
      dydx[i]=b[i].Im;
      dydx[i+NDIM]=-b[i].Re;
    }

  // get back to old matter hamiltonian
  Hmatter.m[0][0].Re-=matterparam;
}

double nuox_propag_in_variable_density(int i,int j)
{
#define equations 2*NDIM
  int a,k,kk;
  int size;
  int nok,nbad;
  double eps=1.e-10;
  double ystart[equations];
  double x1,x2;
  double nstepmin=500.,nstep;
  double hmin=0.0,h;
  double prob;
  double EarthDiam=2.*6371.; // in Km

  complex p,amp;
  matrix Mmass,Mmix,Minducedmass,Mdag,H,Hfin,H1,I;

 // Create mass matrix
  mreset(&Mmass);
  Mmass.m[0][0].Re=NeuOsc.mass2[0];
  Mmass.m[1][1].Re=NeuOsc.mass2[1];
  Mmass.m[2][2].Re=NeuOsc.mass2[2];
 // mprint(Mmass);
  /* Create mixing matrix */
  if(NeuOsc.neutype<0 )
    {
      nuox_set_new_CKM_matrix(NeuOsc.mass2[0],NeuOsc.mass2[1],NeuOsc.mass2[2],NeuOsc.t12,NeuOsc.t13,NeuOsc.t23,-NeuOsc.delta);          // Change delta phase sign for anti-neutrinos
      mcreate(&Mmix,MatOsc.U);      //Changed delta sign!
    }
  else
    mcreate(&Mmix,NeuOsc.U);

  //mprint(Mmass);
  //mprint(Mmix);

  Mdag = mdagger(Mmix);
  H1=mmul(Mmass,Mdag);
//  mprint(H1);
  H=mmul(Mmix,H1);
  mprint(H);

  // Get hamiltonian in matter
  p.Im=0.;
  p.Re=1./(2.*hbarc*1e6*1e-18)/NeuOsc.Energy/1.e9;
  Hmatter=mscale(H,p);
  //mprint(Hmatter);
  // Reset neutrino components
  for (a=0;a<NDIM;a++)
    {
      ystart[a]=ystart[a+NDIM]=0.;
      if(a==i)
	ystart[a]=1;
    }
  
  // Traveled distance
  x1=0.;
  x2=NeuOsc.Length;
 
 // Integration steps from 500 to 4000
  if(NeuOsc.Length <=10.)
    nstep=nstepmin;
  else
    nstep=nstepmin+3500.*NeuOsc.Length/EarthDiam;
  h=(x2-x1)/nstep;

  //Parameters used in odeint
  dxsav=h;
  size =MAXSTEPS*sizeof(double);
  xp=(double*) malloc(size);
  size =equations;
  yp=(double **) malloc((size_t)(size*sizeof(double*)));
  yp[0]=(double *) malloc((size_t)((equations*MAXSTEPS)*sizeof(double)));
  for(kk=1;kk<equations;kk++) yp[kk]=yp[kk-1]+MAXSTEPS;

 // mprint(Hmatter);

  // Integration over traveled path
  odeint(ystart,equations,x1,x2,eps,h,hmin,&nok,&nbad,derivs,rkqs);
  free(xp);
  free(yp[0]);
  free(yp);

  // Finally get the probability!
  prob = SQR(ystart[j])+SQR(ystart[j+NDIM]);
  
  return prob;
}
// ===============================================================================================
// ===============================================================================================
//                              END OF COMPUTATION WITH VARYING DENSITY 
// ===============================================================================================
// ===============================================================================================

/*\docF*/
void nuox_osc_prob_(int *i, int *j, double *proba)
{
*proba = nuox_osc_prob(*i,*j);
}
double nuox_osc_prob(int i, int j)
/* Main entry to calculate single probability between state i to state j */
{
  int k;
  double prob,sin2_32;

  if(NeuOsc.ilevel < 0) /* Non-oscillation case */
      prob = (i==j) ? 1. : 0.;

  if (NeuOsc.ilevel==0) /* Propagation in vacuum */
    {
      prob = nuox_propag(&NeuOsc,i,j);
    }
  
  if (NeuOsc.ilevel==1) /* Propagation in matter (constant density) */
    {
      /* Recompute masses and mixing angles only in case 
         any of the parameters change. Save CPU */
      if(awake_the_master)
	/* Find mass eigenvalues and mixing angles in matter */
	nuox_input_matrix_CKM_matter(i);
      //	nuox_input_matrix_CKM_xing();
      awake_the_master = 0;
      prob = nuox_propag(&MatOsc,i,j);
    }

  if (NeuOsc.ilevel==2) /* Propagation in matter (variable density) */
    {
      prob = nuox_propag_in_variable_density(i,j);
    }
      
  return prob;
}


/***********************************************************************************/
/***********************************************************************************/
/*
                              O B S O L E T E    S T U F F. 
        U S E L E S S   T O   L O O K   B E Y O N D   T H I S   P O I N T
*/
/***********************************************************************************/
/***********************************************************************************/




/**************************************************************************/
/**************************************************************************/
/*    
      VERSION O: USED TO CROSS CHECK, FOR SIMPLE CASES, THE RESULTS 
                 GIVEN BY THE FORMAL THEORY OF NEUTRINO OSCILLATIONS
 */
/**************************************************************************/
/**************************************************************************/

#ifdef CRAP
// The mixing matrix: Parametrization similar to that known as "standard" 
// for the Cabbibo-Kobayashi-Maskawa
void mixing_matrix(double U[3][3]) 
{
  // We assume for the time being that CP is conserved. This simplies things 
  // a good deal

  //First we define the sines and cosines of mixing angles
  c13 = cos(NeuOsc.theta13);
  s13 = sin(NeuOsc.theta13);
  
  c12 = cos(NeuOsc.theta12);
  s12 = sin(NeuOsc.theta12);
  
  c23 = cos(NeuOsc.theta23);
  s23 = sin(NeuOsc.theta23);
  
  //All the elements of the mixing matrix are defined here
  U[0][0]= c12*c13;  /* e1 */
  U[0][1]= s12*c13;  /* e2 */
  U[0][2]= s13;      /* e3 */

  U[1][0]= -s12*c23-c12*s23*s13;  /* mu1 */
  U[1][1]= c12*c23-s12*s23*s13;   /* mu2 */
  U[1][2]= s23*c13;               /* mu3 */

  U[2][0]= s12*s23-c12*c23*s13;  /* tau1 */
  U[2][1]= -c12*s23-s12*c23*s13; /* tau2 */
  U[2][2]= c23*c13;              /* tau3 */
}
#endif

double oscillation_prob(double U[3][3], int i,double L, double E) 
{
  double prob=0.;
  int l,k,iadd;
  double sin2_21,sin2_32,sin2_31;

  // In case we assume CP conservation, the oscillation probability P(l--> k)
  // is the sum of 6 terms. Therefore we have to define the following array
#define MAXTERMS 6
  double term[MAXTERMS];


  switch(i){
    //numu to nue oscillations
  case 12: 
    {  
      l=0;
      k=1;
      break;
    }
    //numu to nutau oscillations
  case 23:
    {
      l=1;
      k=2;
      break;
    }
    //nue to nutau oscillations
  case 13:
    {
      l=0;
      k=2;
      break;
    }
  default: 
    break;
  }

  // We compute the factors containing the oscillatory behaviour
  sin2_21 = sin(1.27*NeuOsc.dm21*L/E)*sin(1.27*NeuOsc.dm21*L/E);
  sin2_32 = sin(1.27*NeuOsc.dm32*L/E)*sin(1.27*NeuOsc.dm32*L/E);
  sin2_31 = sin(1.27*NeuOsc.dm31*L/E)*sin(1.27*NeuOsc.dm31*L/E);

  // We compute the different contributions of each term...
  term[0] = U[l][0]*U[l][0]*U[k][0]*U[k][0];
  term[1] = U[l][1]*U[l][1]*U[k][1]*U[k][1];
  term[2] = U[l][2]*U[l][2]*U[k][2]*U[k][2];

  term[3] = 2.*U[k][0]*U[l][0]*U[k][1]*U[l][1]*(1.-2.*sin2_21);
  term[4] = 2.*U[k][2]*U[l][2]*U[k][1]*U[l][1]*(1.-2.*sin2_32);
  term[5] = 2.*U[k][2]*U[l][2]*U[k][0]*U[l][0]*(1.-2.*sin2_31);

  // Finally sum all the terms to get the oscillation probability
  for(iadd=0;iadd<MAXTERMS;iadd++) 
    prob += term[iadd];

  return prob;
}


#ifdef CRAP
void compute_oscillation_probabilities(double L, double E)
{
  double U[3][3];
  double pnumunue,pnumunutau,pnuenutau;
  
  // Compute mixing matrix 
  mixing_matrix(U);

  // Different oscillation probabilites... (computed in vacuum!!!)

  // numu to nue 
  pnumunue = oscillation_prob(U,12,L,E);
  // numu to nutau 
  pnumunutau = oscillation_prob(U,23,L,E);
  // nue to nutau 
  pnuenutau = oscillation_prob(U,13,L,E);

#ifdef DEBUG
  /* Debugging */
  printf("***************** \n");
  printf("\n");
  printf("Cross check with version 0 \n");
  printf("Nue to Numu %10e \n",pnumunue);
  printf("Nue to Nutau %10e \n",pnuenutau);
  printf("Numu to Nutau %10e \n",pnumunutau);
#endif
}
#endif

#ifdef CRAP
void compute_matter_effects_for_two_neutrino_mixing(double L, double E,int particle)
{
  double prob;
  double x;
  double angle=NeuOsc.theta12;
  double theta_m;
  double factor,lambda;
  double osc_term;
  double sin2thetam;
  double deltam2=NeuOsc.dm21;

  /* Computed for nue to numu oscillations assuming constant density */
  /* particle =1 means neutrinos, particle=-1 means anti-neutrinos */
  x = 0.76e-4 * NeuOsc.density * E/deltam2;
  factor = x-(particle*cos(2*angle));
  
  sin2thetam = sin(2*angle)/sqrt((sin(2*angle)*sin(2*angle))+(factor*factor));
  lambda = L*sqrt(sin(2*angle)*sin(2*angle)+(factor*factor));
  osc_term = 1.27*deltam2*lambda/E;
  prob = sin2thetam*sin2thetam*sin(osc_term)*sin(osc_term); 
#ifdef DEBUG
  printf("  \n");
  printf("Nue to Numu oscillation in matter: %10e \n",prob);
#endif
}
#endif

/**************************************************************************/
/**************************************************************************/
/*    
                    FORMAL THEORY OF NEUTRINO OSCILLATIONS
 */
/**************************************************************************/
/**************************************************************************/

typedef struct M_SPINOR m_spinor;   /* mass eigenstate spinor */
typedef struct W_SPINOR w_spinor;   /* weak eigenstate spinor */

struct M_SPINOR {
  complex component[NDIM];
  double  mass[NDIM];          /* in eV */
};

struct W_SPINOR {
  complex component[NDIM];
};

void create_m_spinor(m_spinor *v, double re_components[NDIM], double masses[NDIM])
{
  int i;
  for(i=0;i<NDIM;i++)
    {
      v->component[i].Re=re_components[i];
      v->component[i].Im=0;
      v->mass[i]=masses[i];         /* eV */
    }
}

void print_m_spinor(char *mess, m_spinor v)
{
  int i;
  printf("m_spinor:  %s      masses: %10e %10e %10e\n",mess,
	 v.mass[0],v.mass[1],v.mass[2]);
  for(i=0;i<NDIM;i++)
    printf("           components:     %i-     %10e  %10e mod(%10e)\n",
	   i,v.component[i].Re,v.component[i].Im,cmod(v.component[i]));
}

void print_w_spinor(char *mess, w_spinor v)
{
  int i;
  printf("w_spinor  %s\n",mess);
  for(i=0;i<NDIM;i++)
    printf("w_spinor   components:     %i-     %10e  %10e mod(%10e)\n",
	   i,v.component[i].Re,v.component[i].Im,cmod(v.component[i]));
}

void print_osc_prob_in_vacuum(double prob[NDIM])
{
  int i;

  printf("   \n");
  printf("***** Oscillation probabilities in vacuum *****\n");

  printf("Probability nu_l to nu_e:    %10e \n",prob[0]);
  printf("Probability nu_l to nu_mu:   %10e \n",prob[1]);
  printf("Probability nu_l to nu_tau:  %10e \n",prob[2]);
}

void compute_osc_prob_in_vacuum(w_spinor w) 
{
  int i;
  double prob[NDIM];
  
  for(i=0;i<NDIM;i++)
    prob[i] = cmod(w.component[i])*cmod(w.component[i]);
#ifdef DEBUG  
  print_osc_prob_in_vacuum(prob);
#endif
}

m_spinor propagate_vacuum(m_spinor v, double p, double L)
/* prograpage neutrino spinor through vacuum over distance L */
/* spinor v, momentum in p (GeV), length L is in meters (km) */
{
  int i;
  m_spinor ve;
  complex expon;
  double Hii;

  ve=v;  
  for(i=0;i<NDIM;i++)
    {
      Hii=v.mass[i]*v.mass[i]/2.0/p/1e9;      /* eV */
      expon.Re=0;
      expon.Im=-Hii*L/(hbarc*1e6*1e-18);
      ve.component[i]=cmul(cexp(expon),ve.component[i]);
    }
  return ve;
}

void create_w_spinor(w_spinor *v, double re_components[NDIM])
{
  int i;
  for(i=0;i<NDIM;i++)
    {
      v->component[i].Re=re_components[i];
      v->component[i].Im=0;
    }
}


void transform_w_m(m_spinor *v, w_spinor w, double masses[NDIM], complex U[NDIM][NDIM])
{
  int i,j;
  complex c,Ustar;
  for(i=0;i<NDIM;i++)
    {
      v->mass[i]=masses[i];         /* eV */
      v->component[i].Re=0;
      v->component[i].Im=0;
      for(j=0;j<NDIM;j++)
	{
	  Ustar=cconj(U[j][i]);
	  c=cmul(w.component[j],Ustar);
	  v->component[i].Re+=c.Re;
	  v->component[i].Im+=c.Im;
	}
    }
}
void transform_m_w(w_spinor *w, m_spinor m,complex U[NDIM][NDIM]) 
{
  int i,j;
  complex c;

  for(i=0;i<NDIM;i++)
    {
      w->component[i].Re=0;
      w->component[i].Im=0;
      for(j=0;j<NDIM;j++)
	{
	  c=cmul(m.component[j],U[i][j]);
	  w->component[i].Re+=c.Re;
	  w->component[i].Im+=c.Im;
	}
    }
}

matrix compute_hamiltonian_vacuum(matrix Mmix,double mass[NDIM],double E,double L) 
{
  
  // In this subroutine we compute: -i*H*L/(hbar *c)
  // where H is defined as (U*H_masses*U^-1)/2*E
  complex c;
  matrix H,H1;
  matrix Mcos,Msin;
  matrix Mmass,MexpH,Mdag;
  matrix Mexp_Hdiag;
  double factor_in_exp=-L/(hbarc*1e6*1e-18);
  
  // Create mass matrix
  mreset(&Mmass);
  Mmass.m[0][0].Re=mass[0]*mass[0];
  Mmass.m[1][1].Re=mass[1]*mass[1];
  Mmass.m[2][2].Re=mass[2]*mass[2];
  
  Mdag = mdagger(Mmix);
  H1=mmul(Mmass,Mdag);
  H=mmul(Mmix,H1);
  
  c.Im=0.;
  c.Re=1./2./NeuOsc.Energy/1.e9;
  H=mscale(H,c);
  
  c.Re=0.;
  c.Im=factor_in_exp; 
  H=mscale(H,c); 
  MexpH=mexp(H);

  return MexpH;
}


void print_probabilities(OscProb P)
{
  printf("  \n");
  printf("***************************************  \n");
  printf("  \n");
  printf("Nue to Nue    --> %10e\n",P.pnuenue);
  printf("Nue to Numu   --> %10e\n",P.pnuenumu);
  printf("Nue to Nutau  --> %10e\n",P.pnuenutau);
  printf("  \n");
  printf("Numu to Nue   --> %10e\n",P.pnumunue);
  printf("Numu to Numu  --> %10e\n",P.pnumunumu);
  printf("Numu to Nutau --> %10e\n",P.pnumunutau);
  printf("  \n");
  printf("Nutau to Nue  --> %10e\n",P.pnutaunue);
  printf("Nutau to Numu  --> %10e\n",P.pnutaunumu);
  printf("Nutau to Nutau --> %10e\n",P.pnutaunutau);

  printf("  \n");
}

OscProb compute_probability_from_amplitude(matrix A) 
{
  matrix P;
  int i,j;
  OscProb Prob;
  
  /* Take the moduli of the amplitudes */
  P=mmod_of_elements(A);
  
  /* Now I fill a new structure that contains all the probabilities for 
     a three family mixing */
  Prob.pnuenue     =mget_prob(P,0,0);
  Prob.pnuenumu    =mget_prob(P,0,1);
  Prob.pnuenutau   =mget_prob(P,0,2);
  
  Prob.pnumunue    =mget_prob(P,1,0);
  Prob.pnumunumu   =mget_prob(P,1,1);
  Prob.pnumunutau  =mget_prob(P,1,2);
  
  Prob.pnutaunue   =mget_prob(P,2,0);
  Prob.pnutaunumu  =mget_prob(P,2,1);
  Prob.pnutaunutau =mget_prob(P,2,2);
#ifdef DEBUG  
  print_probabilities(Prob);
#endif
  return Prob;
} 



void check_probabilities(OscProb *ProbVac)
/* Make sure that the probabilities we computed make some sense. */
{
  double prec=1e-5;
  int ishit=0;

  if(fabs(ProbVac->pnuenue+ProbVac->pnuenumu+ProbVac->pnuenutau-1.0)>prec)
    {
      printf("electron prob not conserved \n");
      ishit++;
    }
  if(fabs(ProbVac->pnumunue+ProbVac->pnumunumu+ProbVac->pnumunutau-1.0)>prec)
    {
      printf("muon prob not conserved %f\n",ProbVac->pnumunue+ProbVac->pnumunumu+ProbVac->pnumunutau);
      ishit++;
    }
  if(fabs(ProbVac->pnutaunue+ProbVac->pnutaunumu+ProbVac->pnutaunutau-1.0)>prec)
    {
      printf("tau prob not conserved %f\n",ProbVac->pnutaunue+ProbVac->pnutaunumu+ProbVac->pnutaunutau);
      ishit++;
    }
  if(ishit)
    {
      printf("Fatal error in physutil.. please contact operator!\n");
      exit(1);
    }
}

OscProb *neutrino_oscillations_vacuum(void)
{
  m_spinor v,vo;
  w_spinor w,wf;
  double comp[NDIM];
  matrix Mmix,Hmat,Hvac;
  int i;
  double DIV=20.;
  OscProb ProbVac,ProbMat;
  OscProb *PVac;

  if(init_parms)
    {
      printf("FATAL: you must first set the masses and mixing of the neutrinos before computing the osc. prob.\n");
      exit(1);
    }

#ifdef OLD_GREAT_STUFF_CANIBALIZED
  comp[0]=0.0;
  comp[1]=1.0;
  comp[2]=0.0;

  create_w_spinor(&w,comp);

  transform_w_m(&v,w,mass,NeuOsc.U);
#ifdef DEBUG
  print_m_spinor("before osc", v);
#endif
  vo=propagate_vacuum(v, NeuOsc.Energy, NeuOsc.Length);
#ifdef DEBUG
  print_m_spinor("after prop",vo);
#endif
  transform_m_w(&wf,vo,U);
#ifdef DEBUG
  print_w_spinor("after osc",wf);
#endif
  /* Finally compute the probabilities */
  compute_osc_prob_in_vacuum(wf);   
#endif

/**********************************************************************************/
/**********************************************************************************/
/* 
                              M A T R I X    F O R M A L I S M 
*/
/**********************************************************************************/
/**********************************************************************************/
  /* Create mixing matrix that relates weak to mass eigenstates */
  mcreate(&Mmix,NeuOsc.U);

/**********************************************************************************/
/**********************************************************************************/
/* 
          N E U T R I N O   P R O P A G A T I O N   I N   V A C U U M  
*/
/**********************************************************************************/
/**********************************************************************************/
#ifdef DEBUG
  printf(" \n");
  printf(" \n");
  printf(" \n");
  printf(" Compute oscillations using the Matrix formalism \n");
#endif
  /* Compute the Hamiltonian in vacuum */
  midentity(&Hvac);
  for(i=0;i<DIV;i++)
    Hvac=mmul(Hvac,compute_hamiltonian_vacuum(Mmix,NeuOsc.mass,NeuOsc.Energy,NeuOsc.Length/DIV));
  ProbVac=compute_probability_from_amplitude(Hvac);
  
  /* check that result makes sense */
  //  check_probabilities(&ProbVac);

  PVac=&ProbVac;

  return PVac;
}
#ifdef LATER
/**********************************************************************************/
/* 
                        F O R T R A N   I N T E R F A C E 
*/
/**********************************************************************************/
/**********************************************************************************/

void input_matrix_ckm_(double *dm32, double *dm21, double *t12,double *t13,double *t23,double *delta)
{
  nuox_input_matrix_CKM(*dm32,*dm21,*t12,*t13,*t23,*delta);
}
  
void input_experimental_parameters__(double *l,double *e,double *rho,int *ntyp)
{
  //  nuox_set_propag_level(3,*rho);
  nuox_set_neutrino(*l,*e,*ntyp);
}

void fneu_oscillations_matter__(int *ineutrino,double *eneutrino,double Proba[3][3])
{
  int i,j;

  OscProb *Prob;
  //  nuox_set_propag_level(3,NeuOsc.density);
  nuox_set_neutrino(NeuOsc.Length,*eneutrino,*ineutrino);
  Prob=neutrino_oscillations_matter();
  
  u.Prob = *Prob;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      Proba[i][j]=u.P[j][i];
}
#endif


void matter_behaviour(int i,double D)
{
  double A,B,C,S,angle;
  double A2,S2,E2,F2;
  double factor,factor1,factor2,factor3;
  double alpha,beta,E,F;
  double DM32,DM21,DM31,t12m,t13m,t23m,delta13m;
  double M1_2,M2_2,M3_2;

  double c13,s13,c132,s132;
  double c12,s12,c122,s122;
  double c23,s23,c232,s232;
  
  double sm12_2,sm13_2,sm23_2;
  double cdeltam,sdeltam;
  double BB;
  double denom;

  int iflag;

  /* I have to compute the hamiltonian eigenvalues and mixing angles in 
     matter 
     Formulae are extracted from H.W. Zaglauer and K.H. Schwarzer 
                                 Z. Phys C 40 (1988) 273   */

  /* Useful values from vacuum mixing angles */
  c13 = cos(NeuOsc.t13);
  s13 = sin(NeuOsc.t13);
  c132 = SQR(c13);
  s132 = SQR(s13);

  c12 = cos(NeuOsc.t12);
  s12 = sin(NeuOsc.t12);
  c122 = SQR(c12);
  s122 = SQR(s12);

  c23 = cos(NeuOsc.t23);
  s23 = sin(NeuOsc.t23);
  c232 = SQR(c23);
  s232 = SQR(s23);

  /*********************************************************************/
  /*                   Effectives masses in matter                     */
  /*********************************************************************/

  /* This is the term: 2* sqrt(2) * G_F * n_e * E  in eV^2 */
  //  D = 7.56e-5*NeuOsc.density*NeuOsc.Energy;
  D *=NeuOsc.dm32;
  if(i<0) /* Anti-neutrinos */
    D = -D;

  A = NeuOsc.dm21 + NeuOsc.dm31 + D;
  A2 = SQR(A);
  
  B = NeuOsc.dm21*NeuOsc.dm31+
    D*(NeuOsc.dm31*c132+NeuOsc.dm21*(c132*c122+s132));

  C = D*NeuOsc.dm21*NeuOsc.dm31*c132*c122;

  angle = (2.*CUBE(A)-9.*A*B+27.*C)/(2.*pow(A2-3.0*B,3./2.));
  angle=acos(angle);
  S = cos((1./3.)*angle);
  S2=SQR(S);

  /* The mass eigenvalues squared are...*/
  factor = NeuOsc.mass2[0] + (A/3.);
  factor1 = sqrt(A2-3*B)*S/3.;
  factor2 = sqrt(A2-3*B)*sqrt(1-S2)/3.;
  M1_2 = factor - factor1-sqrt(3.)*factor2;
  M2_2 = factor - factor1+sqrt(3.)*factor2;
  M3_2 = factor + 2*factor1;

  if(M1_2<0)
    M1_2=0.;

  /* Mass differences in matter */
  DM21=M2_2-M1_2;
  DM31=M3_2-M1_2;
  DM32=M3_2-M2_2;

  /*********************************************************************/
  /* Now we have to define the angles in matter                        */
  /*********************************************************************/

  /* More useful definitions */
  alpha = NeuOsc.mass2[2]*c132+
          NeuOsc.mass2[1]*(c132*c122+s132)+
          NeuOsc.mass2[0]*(c132*s122+s132);

  beta = NeuOsc.mass2[2]*c132*
                  (NeuOsc.mass2[1]*c122+NeuOsc.mass2[0]*s122)+
         NeuOsc.mass2[1]*NeuOsc.mass2[0]*s132;

  E=c13*s13*(NeuOsc.dm31*(M3_2-NeuOsc.mass2[0]-NeuOsc.dm21) - 
             NeuOsc.dm21*(M3_2-NeuOsc.mass2[0]-NeuOsc.dm31)*s122);
  E2 = SQR(E);

  F = (M3_2-NeuOsc.mass2[0]-NeuOsc.dm31)*NeuOsc.dm21*c12*s12*c13;
  F2 = SQR(F);

  /* Now the squared sine of the mixing angles */
  sm12_2 = -(SQR(M2_2)-alpha*M2_2+beta)*DM31;
  sm12_2 = sm12_2/(DM32*(SQR(M1_2)-alpha*M1_2+beta)-
		   DM31*(SQR(M2_2)-alpha*M2_2+beta));

  sm13_2 = (SQR(M3_2)-alpha*M3_2+beta)/DM31/DM32;

  sm23_2 = E2*s232+F2*c232+2*E*F*c23*s23*cos(NeuOsc.delta);
  sm23_2 = sm23_2/(E2+F2);
  //sm23_2=sm23_2/SQR(E+F);
  printf("F,G, sm23= %10.9f %10.9f %10.9f %10.9f\n",F,E,sm23_2,sm12_2);
  /* To avoid stupid crashes due to precision */
  iflag = 0;
  if(sm12_2 > 1.) 
    {
      sm12_2 = 1.;
      iflag =1;
    }
  if(sm13_2 > 1.) 
    {
      sm13_2 = 1.;
      iflag =1;
    }
  if(sm23_2 > 1.) 
    {
      sm23_2 = 1.;
      iflag =1;
    }
  
  if(sm12_2 < 0) 
    {
      sm12_2 = 0;
      iflag =1;
    }
  if(sm13_2 < 0) 
    {
      sm13_2 = 0;
      iflag =1;
    }
  if(sm23_2 < 0) 
    {
      sm23_2 = 0;
      iflag =1;
    }

  if(iflag)printf("fuckup \n");

  t12m = asin(sqrt(sm12_2));
  t13m = asin(sqrt(sm13_2));
  t23m = asin(sqrt(sm23_2));

  
  /*********************************************************************/
  /* And now the CP violation phase in matter !!!!*/
  /*********************************************************************/
  factor3 = 2*E*F*c23*s23*cos(NeuOsc.delta);
  cdeltam = (E2-F2)*cos(NeuOsc.delta)*s23*c23+E*F*(c232-s232);
  sdeltam = -(E2+F2)*sin(NeuOsc.delta)*s23*c23;
  denom = (E2*s232+F2*c232+factor3)*(E2*c232+F2*s232-factor3);
  if(denom<0)
    {
      printf("denom is negative??? %f\n",denom);
    }
  denom=sqrt(denom);

  if(denom>0)
    {
      cdeltam = cdeltam/denom;
      if(cdeltam>1.0)cdeltam=1.0;
      if(cdeltam<-1.0)cdeltam=-1.0;
      if(sdeltam<0)
	delta13m = acos(cdeltam);
      else
	delta13m = - acos(cdeltam);
    }
  else
    delta13m=NeuOsc.delta;


  if(NeuOsc.neutype<0) /* Anti-neutrinos */
    delta13m = - delta13m;
  nuox_set_new_CKM_matrix(M1_2,M2_2,M3_2,t12m,t13m,t23m,delta13m);
}
