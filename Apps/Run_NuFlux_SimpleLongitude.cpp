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

int main(int argc, char * argv[]) {

  // Instantiate DiscreteEarth
  DiscreteEarth d(300.0); // km cell size

  // We set up a given Radiogenic Composition Model for Earth
  // this is done by the DiscreteEarth class
  //  d.SetUniformMantle(d.DepMantle);
  d.SetMantleP1();
  
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
  double dtheta = PIGREEK/180.0;
  double theta_min = 0.0;
  double theta_max = PIGREEK;


  cout << "Number of surface cells along latitude: " << (PIGREEK-0.0)/dtheta << endl;

  ofstream outfile;
  //  outfile.open("Sinogram_uniformmantle.dat");
  outfile.open("Sinogram_pointmantle.dat");

  // Prototype of Radon transform: 
  // - rotate the detector around the patient
  // - collect the 'rays' in a plane for each angle


  // Rotating the Earth
  // Loop over theta (latitude)
  // l = 0-180 degrees
  for(theta = theta_min; theta <= theta_max; theta = theta + dtheta){
    if(int(theta*180.0/PIGREEK) % 10 == 0)
      cout << "Rotation: " << theta*180.0/PIGREEK << endl;


    // Rotate Earth around axis given by the normal to the longitudinal plane
    d.RotateEarth(dtheta, 0, 0, 0);

    // Pick a cell as a reference vector pointing to the local plane of the detector
    Cell_t loc_det = d.GetSurfaceCell(theta, l);

    //    d.PrintCell(loc_det);
    double mag_loc_det = sqrt(loc_det.x*loc_det.x + loc_det.y*loc_det.y + loc_det.z*loc_det.z);// this should be R_e
    // We need to calculate:
    // - Flux vectors with direction parallel to the vector pointing to this cell (collimator approx)
    // - Coordinate of the flux in the local coordinate system of the local plane (inner product with the tangent vector)

    // Loop over all cells along the longitude to collect all 
    // Get list of cells at a given Longitude and theta
    // theta above is taken as center of detector so, we need to go from theta - 90 to theta + 90 degrees
    // But since theta goes from 0 to 180 only, we do this in 2 steps:
    // First 0 to thetamax, then remaining with the 180 degree shifted part
    double theta_l_min = theta - PIGREEK/2.; 
    double theta_l_max = theta + PIGREEK/2.;
    
    for(double theta_l = 0; theta_l <= theta_l_max; theta_l = theta_l + dtheta){
      Cell_t s = d.GetSurfaceCell(theta_l, l);

      // Need to calculate the position of this event in the local (rotated) detector coordinate
      // We create a normalized vector in the local detector plane, which is also in the longitude's plane.

      // 1. We need first a normal vector in the longitude's plane
      // Let vec(P0) = (P0x, P0y, P0z) be a point given in the plane of the longitude
      // let vec(n) = (nx, ny, nz) an orthogonal vector to this plane
      // then vec(P) = (Px, Py, Pz) will be in the plane if (vec(P) - vec(P0)) * vec(n) = 0
      
      // We pick 2 vectors in the plane
      double P0x, P0y, P0z; // given by the longitude
      d.ToCartesian(R_E, 2.0, l, &P0x,&P0y, &P0z);
      double P1x, P1y, P1z;
      d.ToCartesian(R_E, 2.5, l, &P1x,&P1y, &P1z);
      double P2x, P2y, P2z;
      d.ToCartesian(R_E, 3.0, l, &P2x,&P2y, &P2z);

      double v1x = P1x - P0x;
      double v1y = P1y - P0y;
      double v1z = P1z - P0z;
      double v2x = P2x - P0x;
      double v2y = P2y - P0y;
      double v2z = P2z - P0z;
      
      // The cross product will give a vector orthogonal to the plane
      double nx, ny, nz;
      nx = v1y*v2z - v1z*v2y;
      ny = v1z*v2x - v1x*v2z;
      nz = v1x*v2y - v1y*v2x;
      double nMag = sqrt(nx*nx + ny*ny + nz*nz);
      nx /= nMag;
      ny /= nMag;
      nz /= nMag;
      
      // 2. We take the cross product of this normal vector with the vector pointing to the local detector point
      double detpl_tan_x, detpl_tan_y, detpl_tan_z;
      detpl_tan_x = ny*loc_det.z - nz*loc_det.y;
      detpl_tan_y = nz*loc_det.x - nx*loc_det.z;
      detpl_tan_z = nx*loc_det.y - ny*loc_det.x;
      double detpl_tan_mag = sqrt(detpl_tan_x*detpl_tan_x + detpl_tan_y*detpl_tan_y + detpl_tan_z*detpl_tan_z);
      detpl_tan_x /= detpl_tan_mag;
      detpl_tan_y /= detpl_tan_mag;
      detpl_tan_z /= detpl_tan_mag;

      // 3. We take the inner product of the det plane tangent vec with the current surface cell's vector where flux is calculated
      // which should give the local coordinate in the (rotated) detector plane
      double local_coord = detpl_tan_x*s.x + detpl_tan_y*s.y + detpl_tan_z*s.z;


      double FluxU238 = 0; // unit: 
      double meanProb = 0.544; 
      double dx, dy, dz, dr2, mag_d, prod, cos_theta;
      std::vector<Cell_t> cells_l = d.GetCellsLongitude(l);
      for(unsigned int i = 0; i < cells_l.size(); i++){
	if(!d.IsEqual(cells_l[i],s) &&
	   cells_l[i].a238U > 0){
	  // vectorial difference
	  dx = s.x - cells_l[i].x;
	  dy = s.y - cells_l[i].y;
	  dz = s.z - cells_l[i].z;
	  dr2 = dx*dx + dy*dy + dz*dz;
	  mag_d = sqrt(dr2);
	  // Collimator approximation: only accept neutrinos that are parallel
	  // with the local detector's plane normal vector (the vector pointing to the longitude plane)
	  // u*v = |u|*|v|*cos(theta)
	  prod = dx*loc_det.x+dy*loc_det.y+dz*loc_det.z;
	  cos_theta = prod/(mag_loc_det*mag_d);

	  if(fabs(cos_theta - 1.0) < 0.01){
	    FluxU238 += cells_l[i].a238U * cells_l[i].rho / dr2; // [1/kg]*[kg/km3]/[km2] = [1/km5]
	  }
	}
      }
      // apply constants
      //[1/km5] * [km3] * [1] * [1/s] * [1] = [1/(km2*s)]
      FluxU238 *= (d.m_DCell*d.m_DCell*d.m_DCell) * d.U238.Mnu * d.U238.Lamb * meanProb / (4*PIGREEK);
      // cout << theta << "\t" << theta_l*180/PIGREEK << "\t" << local_coord << "\t" << FluxU238/(100000.*100000.*1.e+06) << "\t" << endl;
      outfile << theta << "\t" << local_coord << "\t" << FluxU238/(100000.*100000.*1.e+06) << "\t" << endl;
    }

    for(double theta_l = theta_l_min; theta_l <= 0; theta_l = theta_l + dtheta){
      Cell_t s = d.GetSurfaceCell(theta_l, l);

      // Need to calculate the position of this event in the local (rotated) detector coordinate
      // We create a normalized vector in the local detector plane, which is also in the longitude's plane.

      // 1. We need first a normal vector in the longitude's plane
      // Let vec(P0) = (P0x, P0y, P0z) be a point given in the plane of the longitude
      // let vec(n) = (nx, ny, nz) an orthogonal vector to this plane
      // then vec(P) = (Px, Py, Pz) will be in the plane if (vec(P) - vec(P0)) * vec(n) = 0
      
      // We pick 2 vectors in the plane
      double P0x, P0y, P0z; // given by the longitude
      d.ToCartesian(R_E, 2.0, l, &P0x,&P0y, &P0z);
      double P1x, P1y, P1z;
      d.ToCartesian(R_E, 2.5, l, &P1x,&P1y, &P1z);
      double P2x, P2y, P2z;
      d.ToCartesian(R_E, 3.0, l, &P2x,&P2y, &P2z);

      double v1x = P1x - P0x;
      double v1y = P1y - P0y;
      double v1z = P1z - P0z;
      double v2x = P2x - P0x;
      double v2y = P2y - P0y;
      double v2z = P2z - P0z;
      
      // The cross product will give a vector orthogonal to the plane
      double nx, ny, nz;
      nx = v1y*v2z - v1z*v2y;
      ny = v1z*v2x - v1x*v2z;
      nz = v1x*v2y - v1y*v2x;
      double nMag = sqrt(nx*nx + ny*ny + nz*nz);
      nx /= nMag;
      ny /= nMag;
      nz /= nMag;
      
      // 2. We take the cross product of this normal vector with the vector pointing to the local detector point
      double detpl_tan_x, detpl_tan_y, detpl_tan_z;
      detpl_tan_x = ny*loc_det.z - nz*loc_det.y;
      detpl_tan_y = nz*loc_det.x - nx*loc_det.z;
      detpl_tan_z = nx*loc_det.y - ny*loc_det.x;
      double detpl_tan_mag = sqrt(detpl_tan_x*detpl_tan_x + detpl_tan_y*detpl_tan_y + detpl_tan_z*detpl_tan_z);
      detpl_tan_x /= detpl_tan_mag;
      detpl_tan_y /= detpl_tan_mag;
      detpl_tan_z /= detpl_tan_mag;

      // 3. We take the inner product of the det plane tangent vec with the current surface cell's vector where flux is calculated
      // which should give the local coordinate in the (rotated) detector plane
      double local_coord = detpl_tan_x*s.x + detpl_tan_y*s.y + detpl_tan_z*s.z;


      double FluxU238 = 0; // unit: 
      double meanProb = 0.544; 
      double dx, dy, dz, dr2, mag_d, prod, cos_theta;
      std::vector<Cell_t> cells_l = d.GetCellsLongitude(l);
      for(unsigned int i = 0; i < cells_l.size(); i++){
	if(!d.IsEqual(cells_l[i],s) &&
	   cells_l[i].a238U > 0){
	  // vectorial difference
	  dx = s.x - cells_l[i].x;
	  dy = s.y - cells_l[i].y;
	  dz = s.z - cells_l[i].z;
	  dr2 = dx*dx + dy*dy + dz*dz;
	  mag_d = sqrt(dr2);
	  // Collimator approximation: only accept neutrinos that are parallel
	  // with the local detector's plane normal vector (the vector pointing to the longitude plane)
	  // u*v = |u|*|v|*cos(theta)
	  prod = dx*loc_det.x+dy*loc_det.y+dz*loc_det.z;
	  cos_theta = prod/(mag_loc_det*mag_d);

	  if(fabs(cos_theta - 1.0) < 0.01){
	    FluxU238 += cells_l[i].a238U * cells_l[i].rho / dr2; // [1/kg]*[kg/km3]/[km2] = [1/km5]
	  }
	}
      }
      // apply constants
      //[1/km5] * [km3] * [1] * [1/s] * [1] = [1/(km2*s)]
      FluxU238 *= (d.m_DCell*d.m_DCell*d.m_DCell) * d.U238.Mnu * d.U238.Lamb * meanProb / (4*PIGREEK);
      //      cout << theta << "\t" << theta_l*180/PIGREEK << "\t" << local_coord << "\t" << FluxU238/(100000.*100000.*1.e+06) << "\t" << endl;
      outfile << theta << "\t" << local_coord << "\t" << FluxU238/(100000.*100000.*1.e+06) << "\t" << endl;
    }





  }

  outfile.close();


  // My: 0 --> Glob: +90
  // My: 90 --> Glob: 0
  // My: 180 --> Glob -90

  return 0;
}
