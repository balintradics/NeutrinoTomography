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

// In this example the detector is fixed and the Earth is rotating

int main(int argc, char * argv[]) {

  // Instantiate DiscreteEarth
  DiscreteEarth d(200.0); // km cell size

  // We set up a given Radiogenic Composition Model for Earth
  // this is done by the DiscreteEarth class
  d.SetUniformMantle(d.DepMantle);
  d.SetMantleP1();
  //  d.SetMantleP2();

  // Prepare output file
  ofstream outfile;
  ofstream outfile_det;
  //outfile.open("Sinogram_uniformmantle_Long10_200km.dat");
  //outfile_det.open("Det_sinogram_uniformmantle_Long10_200km.dat");

  outfile.open("Sinogram_pointmantle_Long10_200km.dat");
  outfile_det.open("Det_sinogram_pointmantle_Long10_200km.dat");
  //outfile.open("Sinogram_point2mantle.dat");




  
  //  d.SaveCellsLongitudeToFile(PIGREEK*0.0/180.0);
  //d.SaveSurfaceCellsToFile();
  //d.SaveFluxMap("test_flux_map.dat");

  // We calculate the local Neutrino Flux at each surface position
  // by summing over all contributing cells
  // and evaluating the formula from O. Sramek et al.
  // Flux(r) = (nX * LambX/4pi) * <P> * \int aX(r')*rho(r')/(|r-r'|)^2 dr'
  // units: 1 * [1/s] * 1 * [1/kg] * [kg/m3] * m3 / m2 = [m2 / s]
  // And save it to an output file for each lat, long

  double l = 10*PIGREEK/180.0;
  double l2 = (10+180)*PIGREEK/180.0;
  double theta = 0.0;
  double dtheta = 1.0*PIGREEK/180.0;//d.m_DRad/2.0;
  double theta_min = 0.0;
  double theta_max = 2.0*PIGREEK;

  cout << "Total number of rotations: " << (theta_max - theta_min)/dtheta << endl;

  // Get list of surface cells along longitude - in unrotated coordinates!
  std::vector<Cell_t> surfCells_Long = d.GetSurfaceCellsLongitude(l);



  // Prototype of Radon transform: 
  // - rotate patient and keep the detector fixed 
  // - collect the 'collimated rays' in a plane for each angle

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

  d.CreateDetector(basis1, basis2, normQ, 0);

  // Save ListMode binary file
  std::string outlmfname = std::string("lmfile.dat");
  d.OpenLMF(outlmfname);

  
  // Rotating the Earth
  // Loop over theta (latitude)
  // l = 0-180 degrees
  for(theta = theta_min; theta <= theta_max; theta = theta + dtheta){
    if(int(theta*180.0/PIGREEK) % 10 == 0)
      cout << "Rotation: " << theta*180.0/PIGREEK << endl;

    // Rotate Earth around axis given by the normal to the longitudinal plane
    d.RotateEarth(dtheta, normQ.x, normQ.y, normQ.z);

    // We need to calculate:
    // - Flux vectors with direction parallel to the detector coordinate basis1 vector (collimator approx)
    // - Coordinate of the flux in the detector coordinate system of the longitude's plane (inner product with the basis vector)

    // The detector cells
    std::vector<float> Det1_v(d.m_Det1.size());
 
    // Loop over all cells along the longitude to collect all
    for(int is = 0; is < surfCells_Long.size(); is++){
      // Current cell for which we calculate the flux
      Cell_t s = surfCells_Long[is];
      int detbin = d.GetDetBin(s, basis2);
      // only check half circle
      double rc, thetac, phic;
      d.ToSpherical(s.x, s.y, s.z, &rc, &thetac, &phic);
      //  cout << l << ", " << phic << ", " << d.m_DRad << endl;
      if( fabs(phic-(l+PIGREEK)) > d.m_DRad)
       	continue;
      
      // Need to calculate the position of this point in the fixed detector coordinate system
      double local_coord = s.x * basis2.x + s.y * basis2.y + s.z * basis2.z;

      double FluxU238 = 0; // unit: 
      double meanProb = 0.544; 
      double dx, dy, dz, dr2, mag_d, prod, cos_theta;
      // This gives the already rotated cells!
      std::vector<Cell_t> cells_l = d.GetCellsLongitude(l);
      for(unsigned int i = 0; i < cells_l.size(); i++){
	if(!d.IsEqual(cells_l[i],s)){
	  // vectorial difference
	  dx = s.x - cells_l[i].x;
	  dy = s.y - cells_l[i].y;
	  dz = s.z - cells_l[i].z;
	  dr2 = dx*dx + dy*dy + dz*dz;
	  mag_d = sqrt(dr2);
	  // Collimator approximation: only accept neutrinos that are parallel
	  // with the detector plane basis
	  // u*v = |u|*|v|*cos(theta)
	  prod = dx*basis1.x+dy*basis1.y+dz*basis1.z;
	  cos_theta = prod/(1*mag_d);
	  if(fabs(cos_theta - 1.0) < 0.01){ // 0.05 should correspnd to ~20 degrees resolution...?

	    FluxU238 += cells_l[i].a238U * cells_l[i].rho / dr2; // [1/kg]*[kg/km3]/[km2] = [1/km5]
	  }
	}
      }
      // apply constants
      //[1/km5] * [km3] * [1] * [1/s] * [1] = [1/(km2*s)]
      FluxU238 *= (d.m_DCell*d.m_DCell*d.m_DCell) * d.U238.Mnu * d.U238.Lamb * meanProb / (4*PIGREEK);
      Det1_v[detbin] += FluxU238;
      
      outfile << theta << "\t" << local_coord << "\t" << FluxU238/(100000.*100000.*1.e+06) << "\t" << endl;
      
      //      cout << theta << "\t" << local_coord << "\t" << FluxU238/(100000.*100000.*1.e+06) << "\t" << endl;
    }

    // Save detector cell values as "image"
    for(int idet = 0; idet < Det1_v.size(); idet++){
      outfile_det << Det1_v[idet] << "\t" << flush;
    }
    outfile_det << endl;

  }

  outfile.close();

  outfile_det.close();

  

  d.ProcessLMF(0,0);
  d.CloseLMF();
  
  // My: 0 --> Glob: +90
  // My: 90 --> Glob: 0
  // My: 180 --> Glob -90

  return 0;
}
