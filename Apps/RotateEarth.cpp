#include "stdio.h"
#include "math.h"
#include <iostream>

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

int main(int argc, char * argv[]) {

  // Instantiate DiscreteEarth
  DiscreteEarth d(200.0); // km cell size

  // We set up a given Radiogenic Composition Model for Earth
  // this is done by the DiscreteEarth class
  //d.SetUniformMantle(d.DepMantle);
  d.SetMantleP1();


  d.SaveCellsLongitudeToFile(PIGREEK*10.0/180.0, "Longitude_BeforeRot.dat");

  // Get axis of rotation as a normal of the plane of the Longitude
  Quat4d_t normQ = d.GetNormalToLongitude(PIGREEK*10/180.0);

  // Rotate Earth 
  d.RotateEarth(180.0*PIGREEK/180.0, normQ.x, normQ.y, normQ.z);

  d.SaveCellsLongitudeToFile(PIGREEK*10.0/180.0, "Longitude_AfterRot.dat");


  return 0;
}
