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

  // Given an Earth Cell, its "Depth" depends on from which
  // direction we are looking down at it. This direction is 
  // parametrized by two angles
  
  // To send a neutrino from a given Depth under an angle
  // to a point on the surface, we need to give the
  // - Travel Distance: vectorial difference between the Earth's surface
  // point and the neutrino source
  // - Density profile: along the line of sight for the neutrino
  // This is provided by the DiscreteEarth class
  
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

  // We calculate the Neutrino Oscillation prob for only cells in the plane of the longitude
  // running the oscillation for each to calculate P
  // This is in order to later evaluate the modified formula from O. Sramek et al.:
  // Flux(r) = (nX * LambX/4pi) * \int P(r', r) * aX(r')*rho(r')/(|r-r'|)^2 dr'
  // units: 1 * [1/s] * 1 * [1/kg] * [kg/m3] * m3 / m2 = [m2 / s]
  
  double l = 10*PIGREEK/180.0;
  double l2 = (10+180)*PIGREEK/180.0;
  double theta = 0.0;
  double dtheta = d.m_DRad;
  double theta_min = 0.0;
  double theta_max = PIGREEK;

  // Get list of surface cells along longitude - in unrotated coordinates!
  std::vector<Cell_t> surfCells_Long = d.GetSurfaceCellsLongitude(l);

  cout << "Total cells along longitude: " << surfCells_Long.size() << endl;

  ofstream ofile;
  ofile.open("NuProb_300km_long10deg.txt");

  // Loop over all cells along the longitude to collect all
  for(int is = 0; is < surfCells_Long.size(); is++){
    // Current cell for which we calculate the prob
    Cell_t s = surfCells_Long[is];

    d.PrintCell(s);
    
    // Get list of cells at a given Longitude
    std::vector<Cell_t> cells = d.GetCellsLongitude(l);
    cout << "Looping over " << cells.size() << " longitude cells " << endl;
    double Prob = 1.0; 
    double dx, dy, dz, dr2;
    // Loop over each cell in this slice
    for(unsigned int i = 0; i < cells.size(); i++){
      if(!d.IsEqual(cells[i],s)){

	// Then start from a vector pointing to the original cell
	// and incrementally add a scaled difference vector to it,
	// until it reaches the target
	double Length = d.SetOriginTarget(cells[i], s);
	
	/* neutrino prop in matter */
	nuox_set_propag_level(2,0);
	/* anti-neutrino oscillation */
	nuox_input_matrix_CKM(dm32,dm21,t12,t13,t23,delta);
	nuox_set_neutrino(Length,e,-1);
	Prob=nuox_osc_prob(NU_ELECTRON,NU_ELECTRON);

	ofile << cells[i].x << "\t" << cells[i].y << "\t" << cells[i].z << "\t" << theta << "\t" << l << "\t" << Prob <<  endl;
	//	cout << cells[i].x << "\t" << cells[i].y << "\t" << cells[i].z << "\t" << theta << "\t" << l << "\t" << Prob <<  endl;

      }
    }
  }

  
  // // Loop over theta (latitude + 180 degree to make a full circle)
  // l = l2;
  // for(theta = theta_min; theta <= theta_max; theta = theta + dtheta){
    
  
  //   Cell_t s = d.GetSurfaceCell(theta, l);
  //   d.PrintCell(s);
    
  //   // Get list of cells at a given Longitude
  //   std::vector<Cell_t> cells = d.GetCellsLongitude(l);
    
  //   double Prob = 1.0; 
  //   double dx, dy, dz, dr2;
  //   // Loop over each cell in this slice
  //   for(unsigned int i = 0; i < cells.size(); i++){
  //     if(!d.IsEqual(cells[i],s) &&
  // 	 cells[i].a238U > 0){
  // 	// Then start from a vector pointing to the original cell
  // 	// and incrementally add a scaled difference vector to it,
  // 	// until it reaches the target
  // 	double Length = d.SetOriginTarget(cells[i], s);
	
  // 	/* neutrino prop in matter */
  // 	nuox_set_propag_level(2,0);
  // 	/* anti-neutrino oscillation */
  // 	nuox_input_matrix_CKM(dm32,dm21,t12,t13,t23,delta);
  // 	nuox_set_neutrino(Length,e,-1);
  // 	Prob=nuox_osc_prob(NU_ELECTRON,NU_ELECTRON);

  // 	ofile << cells[i].x << "\t" << cells[i].y << "\t" << cells[i].z << "\t" << theta << "\t" << l << "\t" << Prob <<  endl;

  //     }
  //   }
  // }





  ofile.close();

  return 0;
}
