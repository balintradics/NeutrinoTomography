#include "stdio.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <vector>

#include "../EarthModel/DiscreteEarth.h"
#include "../NeutrinoOsc/neutrino_osc.h"

#define PIGREEK 3.141592654

// Need to re-declare static variables
// that are used in the classes
double DiscreteEarth::m_Dx;
double DiscreteEarth::m_Dy;
double DiscreteEarth::m_Dz;
double DiscreteEarth::m_PathLength;
Cell_t DiscreteEarth::m_Ocell;

using namespace std;

int main(int argc, char * argv[]) {

  // double Angles (from Koike and Sato!)
  double t12=33.0/180.0*PIGREEK;
  double t13=asin(sqrt(0.025));
  double t23 = PIGREEK/4.0;
  // deltaCP
  double delta=-90./180.0*PIGREEK;
  // Differences in mass squared
  double dm32=2.32e-3;
  double dm21=7.59e-5;
  double e = 0.003; // neutrino energy GeV
  
  // Instantiate DiscreteEarth
  DiscreteEarth d(300.0); // km cell size

  // Optionally, we set up a given Radiogenic Composition Model for Earth
  // this is done by the DiscreteEarth class
  d.SetUniformMantle(d.DepMantle);

  //ofstream ofile;
  // ofile.open("AvCellProb.txt");

  // longitude
  double l = 10*PIGREEK/180.0;
  double l2 = (10+180)*PIGREEK/180.0;
  // latitude
  double theta = 45.0*PIGREEK/180.0;

  Cell_t s = d.GetSurfaceCell(theta, l);
  cout << "Surface: " << s.x << ", " << s.y << ", " << s.z << endl;
  double xo, yo, zo;
  d.ToCartesian(R_E-1000, theta, l2, &xo, &yo, &zo);
  Cell_t o = d.GetCell(xo, yo, zo);
  cout << "Origin: " << o.x << ", " << o.y << ", " << o.z << endl;

  double Prob = 1.0; 
  double meanProb = 0.0;
  int N = 1000;
  for(int i = 0; i < N;i++){
    Cell_t oc; // the random cell position
    oc.x = o.x + d.GetRandomNumber() * d.m_DCell;
    oc.y = o.y + d.GetRandomNumber() * d.m_DCell;
    oc.z = o.z + d.GetRandomNumber() * d.m_DCell;
    double Length = d.SetOriginTarget(oc, s);
    //    cout << oc.x << ", " << oc.y << ", " << oc.z << " : " << Length << endl;
    /* neutrino prop in matter */
    nuox_set_propag_level(2,0);
    /* anti-neutrino oscillation */
    nuox_input_matrix_CKM(dm32,dm21,t12,t13,t23,delta);
    nuox_set_neutrino(Length,e,-1);
    Prob=nuox_osc_prob(NU_ELECTRON,NU_ELECTRON);
    meanProb += Prob;
  }
  meanProb /= (double)N;
    
  cout << meanProb << endl;

  return 0;
}

