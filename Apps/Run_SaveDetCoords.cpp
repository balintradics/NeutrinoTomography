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
  d.SetUniformMantle(d.DepMantle);
  //d.SetMantleP1();
  //  d.SetMantleP2();

  // Prepare output file
  ofstream outfile;
  outfile.open("DetectorCoords.dat");


  double l = 10*PIGREEK/180.0;


  // Get list of surface cells along longitude - in unrotated coordinates!
  std::vector<Cell_t> surfCells_Long = d.GetSurfaceCellsLongitude(l);

  // Get axis of rotation as a normal vector of the plane of the Longitude
  Quat4d_t normQ = d.GetNormalToLongitude(l);

  // Get detector plane coordinate system's basis vectors
  Quat4d_t basis1;// X - axis rotated into the Longitude plane
  Quat4d_t basis2;// Z - axis (0, 0, 1)
  d.GetLongitudePlaneBasis(l, &basis1, &basis2);

  cout << "Detector Basis1: " << endl;
  d.PrintQ(basis1);
  cout << "Detector Basis2: " << endl;
  d.PrintQ(basis2);

  d.CreateDetector(basis1, basis2);
  //  d.PrintDetector();

  for(int i = 0; i < d.m_Det1.size(); i++){
  
    outfile <<  d.m_Det1[i].x << "\t" << d.m_Det1[i].y << "\t" << d.m_Det1[i].z <<  endl;
  }
  
  // Loop over all cells along the longitude to collect all
  for(int is = 0; is < surfCells_Long.size(); is++){
    // Current cell for which we calculate the flux
    Cell_t s = surfCells_Long[is];

    
    // only check half circle
    double rc, thetac, phic;
    d.ToSpherical(s.x, s.y, s.z, &rc, &thetac, &phic);
    if( fabs(phic-(l+PIGREEK)) > d.m_DRad)
       	continue;

    
    
    // Need to calculate the position of this point in the fixed detector coordinate system
    // Basically, inner product of the cell coordinate with the z - axis in 3D
    //    double z_coord = s.x * basis2.x + s.y * basis2.y + s.z * basis2.z;
    
    Quat4d_t qdet = d.GetDetCoord(s, basis2);
    d.PrintQ(qdet);
    

  }

  outfile.close();
  
  return 0;
}
