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
  DiscreteEarth d(200.0); // km cell size

  // We set up a given Radiogenic Composition Model for Earth
  // this is done by the DiscreteEarth class
  d.SetUniformMantle(d.DepMantle);
  //d.SetMantleP1();
  //  d.SetMantleP2();

  // Prepare output file
  ofstream outfile;
  ofstream outfile2; // fake detector coords
  
  outfile.open("DetectorCoords_10deg_200km.dat");
  outfile2.open("DetectorCoords_10deg_200km_fake.dat");

  double l = 10*PIGREEK/180.0;
  double l2 = (10+180)*PIGREEK/180.0;
  double theta = 0.0;
  double dtheta = 15.0*PIGREEK/180.0;//d.m_DRad/2.0;
  double theta_min = 0.0;
  double theta_max = 2.0*PIGREEK;
  cout << "Total number of rotations: " << (theta_max - theta_min)/dtheta << endl;

  
  // Get list of surface cells along longitude - in unrotated coordinates!
  std::vector<Cell_t> surfCells_Long = d.GetSurfaceCellsLongitude(l);

  // Get axis of rotation as a normal vector of the plane of the Longitude
  Quat4d_t normQ = d.GetNormalToLongitude(l);

  // Get detector plane coordinate system's basis vectors
  Quat4d_t basis1;// X - axis rotated into the Longitude plane
  Quat4d_t basis2;// Z - axis (0, 0, 1)
  d.GetLongitudePlaneBasis(l, &basis1, &basis2);

  int GlobalDetID1 = 0;
  int GlobalDetID2 = 1000000;// fake detectors opposite shifted by 1e+06
  
  for(theta = theta_min; theta <= theta_max; theta = theta + dtheta){
    if(int(theta*180.0/PIGREEK) % 10 == 0)
      cout << "Rotation: " << theta*180.0/PIGREEK << endl;

    Quat4d_t rbasis1 = d.RotateQuaternion(basis1, normQ, theta);
    Quat4d_t rbasis2 = d.RotateQuaternion(basis2, normQ, theta);
    
    cout << "Detector Basis1: " << endl;
    d.PrintQ(rbasis1);
    cout << "Detector Basis2: " << endl;
    d.PrintQ(rbasis2);
    cout << "--------------------------" << endl;
  
    d.CreateDetector(rbasis1, rbasis2, normQ, theta);
    //  d.PrintDetector();
    
    // SAFIR form for crystal map
    // -----------------------
    // * ring = ringNumber (usually along z-direction)
    // * crystal = crystalNumber (usually around phi direction)
    // * layer = layerNumber, not used at the moment, keep it 0
    // * x, y, z = detector coordinate in mm
    // * angleRad = orientation in phi direction of the normal of the crystal surface, pointing outwardsin
    // * crystalLengthMM = length of the crystal, used to emulate the DOI effect on the system matrix (only important if numberOfPointsPerCrystalDOI is set to > 1)
    
    //#ring   #detector       #layer  x       y       z       angleRad        crystalLengthMM detectorPosition
    //0       0       0       1.570000000000000e+01   6.405000000000000e+01   -1.680000000000000e+01  0.000000000000000e+00   1.300000000000000e+01   0
    //0       1       0       1.350000000000000e+01   6.405000000000000e+01   -1.680000000000000e+01  0.000000000000000e+00   1.300000000000000e+01   0
    //-----------------------
    
    for(int i = 0; i < d.m_Det1.size(); i++){  
      outfile << 0 << "\t" << GlobalDetID1 << "\t" << 0 << "\t" << d.m_Det1[i].x*1e+06 << "\t" << d.m_Det1[i].y*1e+06 << "\t" << d.m_Det1[i].z*1e+06 <<  "\t" << l << "\t" << 1 << "\t" << 0 << endl;
      //    outfile <<  d.m_Det2[i].x << "\t" << d.m_Det2[i].y << "\t" << d.m_Det2[i].z <<  endl;
      GlobalDetID1++;
    }

    for(int i = 0; i < d.m_Det2.size(); i++){  
      outfile2 << 0 << "\t" << GlobalDetID2 << "\t" << 0 << "\t" << d.m_Det2[i].x*1e+06 << "\t" << d.m_Det2[i].y*1e+06 << "\t" << d.m_Det2[i].z*1e+06 <<  "\t" << l << "\t" << 1 << "\t" << 0 << endl;
      //    outfile <<  d.m_Det2[i].x << "\t" << d.m_Det2[i].y << "\t" << d.m_Det2[i].z <<  endl;
      GlobalDetID2++;
    }
    
    // also save the "fake" mirror detector
    // for(int i = 0; i < d.m_Det2.size(); i++){
    //   outfile << 0 << "\t" << d.m_Det1.size()+i << "\t" << 0 << "\t" << d.m_Det2[i].x*1e+06 << "\t" << d.m_Det2[i].y*1e+06 << "\t" << d.m_Det2[i].z*1e+06 <<  "\t" << l2 << "\t" << 1 << "\t" << 0 << endl;
    //   //    outfile <<  d.m_Det2[i].x << "\t" << d.m_Det2[i].y << "\t" << d.m_Det2[i].z <<  endl;
    // }
    
  }


  outfile.close();
  outfile2.close();
  
  return 0;
}
