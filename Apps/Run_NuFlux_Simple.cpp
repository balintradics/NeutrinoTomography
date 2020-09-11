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
  DiscreteEarth d(500.0); // km cell size

  // We set up a given Radiogenic Composition Model for Earth
  // this is done by the DiscreteEarth class
  d.SetUniformMantle(d.DepMantle);

  //  d.SaveActivityMap2DToFile();
  //  d.SaveCellsLongitudeToFile(PIGREEK*45.0/180.0);
  //  d.SaveSurfaceCellsToFile();
  d.SaveFluxMap();


  // We calculate the Neutrino Flux by summing over all cells
  // and simply evaluate the formula from O. Sramek et al.
  // Flux(r) = (nX * LambX/4pi) * <P> * \int aX(r')*rho(r')/(|r-r'|)^2 dr'
  // units: 1 * [1/s] * 1 * [1/kg] * [kg/m3] * m3 / m2 = [m2 / s]
  
  // Get a surface cell at (theta phi)
  Cell_t s = d.GetSurfaceCell(0.0, PIGREEK/2);
  d.PrintCell(s);


  // Integrate contributions from all other cells
  double FluxU238 = 0; // unit: 
  double meanProb = 0.544; 
  double dx, dy, dz, dr2;
  for(unsigned int i = 0; i < d.m_NCells; i++){
    if(!d.IsEqual(d.m_EarthCells[i],s)){
      // vectorial difference
      dx = s.x - d.m_EarthCells[i].x;
      dy = s.y - d.m_EarthCells[i].y;
      dz = s.z - d.m_EarthCells[i].z;
      dr2 = dx*dx + dy*dy + dz*dz;

      FluxU238 += d.m_EarthCells[i].a238U * d.m_EarthCells[i].rho / dr2; // [1/kg]*[kg/km3]/[km2] = [1/km5]
    }
  }
  // apply constants
  //[1/km5] * [km3] * [1] * [1/s] * [1] = [1/(km2*s)]
  FluxU238 *= (d.m_DCell*d.m_DCell*d.m_DCell) * d.U238.Mnu * d.U238.Lamb * meanProb / (4*PIGREEK);

  std::cout << "Flux U238 : " << FluxU238 << " aneutrinos / (km2*s)" << " = " << FluxU238/(100000.*100000.*1.e+06) << " aneutrinos / (cm2*us) " <<  std::endl;

  return 0;
}
