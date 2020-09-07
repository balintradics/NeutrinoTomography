#include "stdio.h"
#include "math.h"

#include "neutrino_osc.h"
//#include "dstconst.h"
//#include "evstore.h"

//#define DMAX 732
//#define EMAX 7.5
#define SHOBAS

#ifdef MEDBAS
#define DMAX 2900
#define EMAX 30
#endif

#ifdef SHOBAS
#define DMAX 732
#define EMAX 7.5
#endif

#define PIGREEK 3.141592654

void plot_osc_vs_l(double rho, char *fname, int inu, int jnu)
{
  // double Angles (from Koike and Sato!)
  //  double t12=22.5/180.0*PIGREEK;
  double t12=33.0/180.0*PIGREEK;

  double t13=asin(sqrt(0.025));
  double t23 = PIGREEK/4.0;
  //  double delta=33.0/180.0*PIGREEK;
  double delta=-90./180.0*PIGREEK;

  // Differences in mass squared
  double dm32=2.32e-3;
  double dm21=7.59e-5;


  double l,e,acp[2],at[2];
  double pnuenumu[2],pnumunue[2],panueanumu[2],panumuanue[2];
  int imatter;
  int ie;

  FILE *f;
  f=fopen(fname,"w");
  double R_E=6371.; // in Km
  double R_OuterCore = 3480.;
// travelling from the Earth outer core to the surface
  for(l=R_E-R_OuterCore; l>0.1;l-=100.0)
  //  for(l=R_E-R_OuterCore; l>=R_E-R_OuterCore;l-=1.0)
    {
       printf("l: %10.2f km\n",l);

      e=0.03; // GeV
      fprintf(f,"%10.2f ",l);
      for(imatter=0;imatter<2;imatter++)
	{
	  if(imatter==0)
	    {
	      /* in vacuuo */
	      nuox_set_propag_level(0,0.);
	    }
	  else
	    {
	      /* in matter */
	      nuox_set_propag_level(2,rho);
	    }
	  
	  /* neutrino oscillation */
	  nuox_input_matrix_CKM(dm32,dm21,t12,t13,t23,delta);
	  //      nuox_dump_osc_matrix();
	  nuox_set_neutrino(l,e,1);
//	 nuox_set_neutrino_costheta(0.);

	  pnuenumu[imatter]=nuox_osc_prob(inu,jnu);
	 // pnumunue[imatter]=nuox_osc_prob(jnu,inu);
	  /* anti-neutrino oscillation */
	  nuox_input_matrix_CKM(dm32,dm21,t12,t13,t23,delta);
	  nuox_set_neutrino(l,e,-1);
	 // nuox_set_neutrino_costheta(0.999);

	  panueanumu[imatter]=nuox_osc_prob(inu,jnu);
	  panumuanue[imatter]=nuox_osc_prob(jnu,inu);
	  
	  fprintf(f,"%10.4e  %10.4e %10.4e %10.4e      ",
		  pnuenumu[imatter],pnumunue[imatter],
		  panueanumu[imatter],panumuanue[imatter]);	  
	}
      fprintf(f,"\n");
    }
  fclose(f);
}


void plot_cp_violation(double l, double rho, char *fname, int inu, int jnu, double t13, double delta)
{
#define SOL3
#ifdef SOL1
  double t12=22.5/180.0*PIGREEK;
  double t13=asin(sqrt(0.025));
  double t23 = PIGREEK/4.0;
  double delta=0./180.0*PIGREEK;
  double dm32=3e-3;
  double dm21=1e-6;
#endif
#ifdef SOL2
  double t12=22.5/180.0*PIGREEK;
  double t13=asin(sqrt(0.025));
  double t23 = PIGREEK/4.0;
  double delta=0./180.0*PIGREEK;
  double dm32=3e-3;
  double dm21=1e-4;
#endif
#ifdef SOL3
  double t12=PIGREEK/4.0;
  //  double t13=asin(sqrt(0.05))/2.;
  double t23 = PIGREEK/4.0;
  //  double delta=-90./180.0*PIGREEK;
  double dm32=3.e-3;
  double dm21=1.e-04;//7.e-5;
#endif

  double e,acp[2],at[2];
  double pnuenumu[2],pnumunue[2],panueanumu[2],panumuanue[2];
  int imatter;
  int ie;

  FILE *f;

  nuox_set_propag_level(0,0.);
  nuox_input_matrix_CKM(dm32,dm21,t12,t13,t23,0.);
  nuox_set_neutrino(700.,10.,1);
  e=nuox_osc_prob(inu,jnu);
  printf("prob %10.4e\n",e);
  nuox_input_matrix_CKM(dm32,dm21,t12,t13,t23,10.);
  nuox_set_neutrino(700.,10.,1);
  e=nuox_osc_prob(inu,jnu);
  printf("prob %10.4e\n",e);

  f=fopen(fname,"w");

  for(ie=0;ie<6400.;ie++)
    {
      e=pow(10.,(-1.0+0.000625*ie));
      fprintf(f,"%10.2f ",e);
      for(imatter=0;imatter<2;imatter++)
	{
	  if(imatter==0)
	    {
	      /* in vacuuo */
	      nuox_set_propag_level(0,0.);
	    }
	  else
	    {
	      /* in matter */
	      nuox_set_propag_level(2,rho);
	    }

	  /* neutrino oscillation */
	  nuox_input_matrix_CKM(dm32,dm21,t12,t13,t23,delta);
	  //      nuox_dump_osc_matrix();
	  nuox_set_neutrino(l,e,1);
	  nuox_set_neutrino_costheta(PIGREEK*90.0/180.0);
	  pnuenumu[imatter]=nuox_osc_prob(inu,jnu);
	  pnumunue[imatter]=nuox_osc_prob(jnu,inu);
	  //	  if(imatter)
	    //	    printf("dm  %10.4e  %10.5e  %10.5e \n",e,MatOsc.dm21,MatOsc.dm32);
	  /* anti-neutrino oscillation */
	  nuox_input_matrix_CKM(dm32,dm21,t12,t13,t23,delta);
	  nuox_set_neutrino(l,e,-1);
	  nuox_set_neutrino_costheta(PIGREEK*90.0/180.0);
	  panueanumu[imatter]=nuox_osc_prob(inu,jnu);
	  panumuanue[imatter]=nuox_osc_prob(jnu,inu);
	  
	  /* CP assymmetry */
	  acp[imatter]=(pnuenumu[imatter]-panueanumu[imatter]);
	  // (pnuenumu[imatter]+panueanumu[imatter]);
	  
	  /* T violation */
	  at[imatter]=(pnuenumu[imatter]-pnumunue[imatter]);
	    ///(pnuenumu[imatter]+pnumunue[imatter]);

	  fprintf(f,"%10.4e  %10.4e %10.4e %10.4e %10.5e %10.5e     ",
		  pnuenumu[imatter],pnumunue[imatter],
		  panueanumu[imatter],panumuanue[imatter],
		  acp[imatter],at[imatter]);

	}
      fprintf(f,"\n");
    }
  fclose(f);
}



int main(int argc, char * argv[]) {


#ifdef SHOBAS
#define DIST 732
#define EMU 7.5
#define DEN 2.8
#endif

#ifdef MEDBAS
#define DIST 2900
#define EMU 30
#define DEN 3.2
#endif

  OscProb *ProbVac;
  OscProb *ProbMat;


double prob,probme,probem,probame,probaem;

  // double Angles 
  //  double t12=3.141592/2.0;
  //  double t13=0.15878;
  double t12;
  double t13;
  double t23;
  double delta=0;//+PIGREEK/2.0;

  // Differences in mass squared
  double dm32=3.e-3;
  double dm21=1.e-4;

  // Experimental parameters
  double l=DIST;
  double e=5.;
  double d=DEN;
  double neve,lresc;
  int ntyp;
  double diff,le;

 // plot_cp_violation(730.,2.8,(char*)("osc730_a.data"),NU_ELECTRON,NU_MUON,asin(sqrt(0.015)),90./180.0*PIGREEK);
 // plot_cp_violation(730.,2.8,(char*)("osc730_d.data"),NU_ELECTRON,NU_MUON,asin(sqrt(0.015)),-90./180.0*PIGREEK);

  plot_osc_vs_l(2.8, (char*)("losc_test.data"), NU_ELECTRON,NU_ELECTRON);


return 0;


}
