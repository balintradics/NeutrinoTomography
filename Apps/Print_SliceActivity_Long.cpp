#include "stdio.h"
#include "math.h"
#include <iostream>
#include <fstream>

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

// Calculating the Sinogram from a Slice of the Earth along a longitude
// using a mean oscillation probability of 0.544

int main(int argc, char * argv[]) {

  // Instantiate DiscreteEarth
  DiscreteEarth d(100.0); // km cell size

  // We set up a given Radiogenic Composition Model for Earth
  // this is done by the DiscreteEarth class
  //d.SetUniformMantle(d.DepMantle);
  //  d.SetMantleP1();
  d.SetMantleP2();


  d.SaveCellsLongitudeToFile(PIGREEK*10.0/180.0, "Activity_MantleP2_long10_100km.dat");

  return 0;

}
