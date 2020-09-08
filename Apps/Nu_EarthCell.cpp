#include "stdio.h"
#include "math.h"
#include <iostream>

#include "../EarthModel/DiscreteEarth.h"
#include "../NeutrinoOsc/neutrino_osc.h"

#define PIGREEK 3.141592654

// Need to re-declare static variables
// that are used in the classes
float DiscreteEarth::m_Dx;
float DiscreteEarth::m_Dy;
float DiscreteEarth::m_Dz;
float DiscreteEarth::m_PathLength;
Cell_t DiscreteEarth::m_Ocell;

int main(int argc, char * argv[]) {


// Given an Earth Cell, its "Depth" depends on from which
// direction we are looking down at it. This direction is 
// parametrized by an angle (later 2 angles for 3D).

// To send a neutrino from a given Depth under an angle
// to a point on the surface, we need to give the
// - Travel Distance: vectorial difference between the Earth's surface
// point and the neutrino source
// - Density profile: along the line of sight for the neutrino

// double Angles (from Koike and Sato!)
double t12=33.0/180.0*PIGREEK;
double t13=asin(sqrt(0.025));
double t23 = PIGREEK/4.0;
// deltaCP
double delta=-90./180.0*PIGREEK;
// Differences in mass squared
double dm32=2.32e-3;
double dm21=7.59e-5;

 
 DiscreteEarth d(100.0); // km cell size

 // Get a random cell
 //Cell_t c = d.GetRandomCell();
 //  Cell_t c = d.GetCell(0.0, 0.0, 0.0);
 Cell_t c = d.GetCell(-2171, 0, -4000);
 d.PrintCell(c);

 // Get a surface cell at (theta phi)
 Cell_t s = d.GetSurfaceCell(0.0, PIGREEK/2);
 d.PrintCell(s);


 // Then start from a vector pointing to the original cell
 // and incrementally add a scaled difference vector to it,
 // until it reaches the target
 float Length = d.SetOriginTarget(c, s);


 // double R_OuterCore = 3480.;


  // Anti-Neutrino energy
  double e = 0.03; // GeV

  //  double depth = R_E - R_LM;

  double prob;

  /* neutrino prop in matter */
  nuox_set_propag_level(2,0);
  /* anti-neutrino oscillation */
  nuox_input_matrix_CKM(dm32,dm21,t12,t13,t23,delta);
  nuox_set_neutrino(Length,e,-1);
  prob=nuox_osc_prob(NU_ELECTRON,NU_ELECTRON);

  std::cout << "Probability: " << prob << std::endl;
  return 0;
}
