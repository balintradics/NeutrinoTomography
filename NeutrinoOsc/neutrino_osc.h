#include "thematrix.h"

#define NDIM 3
#define NU_ELECTRON 0
#define NU_MUON     1
#define NU_TAU      2

// This is the final ouput: the probabilities for all the oscillations processes
typedef struct OSCPROB OscProb;
struct OSCPROB {
  double pnuenue;
  double pnuenumu;
  double pnuenutau;

  double pnumunue;
  double pnumunumu;
  double pnumunutau;
  
  double pnutaunue;
  double pnutaunumu;
  double pnutaunutau;
};

union{
OscProb Prob;
double P[3][3];
}u;

// In this structure we keep all the relevant input information
struct NEUOSC_REC {

  int ilevel;

//Mass differences squared
  double dm32;
  double dm21;
  double dm31;
  double mass[3];
  double mass2[3];    /* mass squared */

  //Angles 
  double t12;
  double t13;
  double t23;
  double delta;

  double sin22t12, sin22t13, sin22t23;
  double sin2t12, sin2t13, sin2t23;

// Length and energy 
  double Length; 
  double Energy;
  
  // Cosine theta for atmospheric neutrinos
  double CosTheta;

// Medium density 
  double density;

//Neutrino type
  int neutype; // 1=neutrinos, -1 means anti-neutrinos

  complex U[NDIM][NDIM];
  complex Ukla[NDIM][NDIM][NDIM];

};

typedef struct NEUOSC_REC NeutrinoOsc;

extern NeutrinoOsc NeuOsc;
extern NeutrinoOsc MatOsc;

// Word needed to avoid computing twice the same hamiltonian in matter
extern int Recompute_Hamiltonian;
void nuox_set_hamiltonian(int i);
void nuox_set_neutrino(double l,double e,int ntyp);
void nuox_input_matrix_CKM(double dm32,double dm21, double t12,double t13,double t23,double delta);
void nux_CKM_from_U(NeutrinoOsc *a);
void nuox_dump_mix_matrix(char *mess, complex U[NDIM][NDIM]);
void nuox_dump_osc_matrix();
void nuox_set_propag_level(int ilevel, double rho);
double nuox_osc_prob(int i, int j);

void compute_oscillation_probabilities(double L, double E); 
OscProb *neutrino_oscillations_vacuum(void);
OscProb *neutrino_oscillations_matter(void);


void nuox_set_neutrino_(double *l,double *e,int *ntyp);
void nuox_set_neutrino_costheta( double );
void nuox_input_matrix_ckm_(double *dm32,double *dm21, double *t12,double *t13,double *t23,double *delta);
void nuox_set_propag_level_(int *ilevel, double *rho);
void nuox_osc_prob_(int *i, int *j, double *);


