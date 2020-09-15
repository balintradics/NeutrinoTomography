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
  
  d.SaveCellsLongitudeToFile(PIGREEK*0.0/180.0);
  //d.SaveSurfaceCellsToFile();
  //d.SaveFluxMap("test_flux_map.dat");

  // We calculate the local Neutrino Flux at each surface position
  // by summing over all contributing cells
  // and evaluating the formula from O. Sramek et al.
  // Flux(r) = (nX * LambX/4pi) * <P> * \int aX(r')*rho(r')/(|r-r'|)^2 dr'
  // units: 1 * [1/s] * 1 * [1/kg] * [kg/m3] * m3 / m2 = [m2 / s]
  // And save it to an output file for each lat, long

  double r, theta, phi;
  // d.ToSpherical(0, 0, 6700, &r, &theta, &phi );
  // cout << "0, 0, 6700: " << r << "\t" << theta << "\t" << phi << endl;

  ofstream outfile;
  outfile.open("Flux_map.dat");
  double TotFluxU238 = 0;
  // Loop over all Surface cell
  for(int i = 0; i < d.m_NSurfCells;i++){
    Cell_t cell = d.m_SurfCells[i];
    
    //d.PrintCell(cell);

    d.ToSpherical(cell.x, cell.y, cell.z, &r, &theta, &phi );
    
    // Integrate contributions from all other cells
    double FluxU238 = 0; 
    double meanProb = 0.544; 
    double dx, dy, dz, dr2;
    for(unsigned int i = 0; i < d.m_NCells; i++){
      if(!d.IsEqual(d.m_EarthCells[i],cell)){
	// vectorial difference
	dx = cell.x - d.m_EarthCells[i].x;
	dy = cell.y - d.m_EarthCells[i].y;
	dz = cell.z - d.m_EarthCells[i].z;
	dr2 = dx*dx + dy*dy + dz*dz;
	
	FluxU238 += d.m_EarthCells[i].a238U * d.m_EarthCells[i].rho / dr2; // [1/kg]*[kg/km3]/[km2] = [1/km5]
      }
    }
    // apply constants
    //[1/km5] * [km3] * [1] * [1/s] * [1] = [1/(km2*s)]
    FluxU238 *= (d.m_DCell*d.m_DCell*d.m_DCell) * d.U238.Mnu * d.U238.Lamb * meanProb / (4*PIGREEK);
    TotFluxU238 += FluxU238;
    outfile << -theta + PI/2 << "\t" << phi << "\t" << FluxU238/(100000.*100000.*1.e+06) << endl;

  }
  std::cout << "Mean Flux U238 : " << TotFluxU238/d.m_NSurfCells << " aneutrinos / (km2*s)" << " = " << TotFluxU238/(100000.*100000.*1.e+06)/d.m_NSurfCells << " aneutrinos / (cm2*us) " <<  std::endl;

  // My: 0 --> Glob: +90
  // My: 90 --> Glob: 0
  // My: 180 --> Glob -90

  return 0;
}
