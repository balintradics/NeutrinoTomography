#include "stdio.h"
#include "math.h"
#include <iostream>
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
  
  
  // Instantiate DiscreteEarth
  DiscreteEarth d(100.0); // km cell size

  // We set up a given Radiogenic Composition Model for Earth
  // this is done by the DiscreteEarth class
  d.SetUniformMantle(d.DepMantle);

  // We calculate the Neutrino Flux by summing over all cells
  // running the oscillation for each to calculate P
  // evaluate the modified formula from O. Sramek et al.:
  // Flux(r) = (nX * LambX/4pi) * \int P(r', r) * aX(r')*rho(r')/(|r-r'|)^2 dr'
  // units: 1 * [1/s] * 1 * [1/kg] * [kg/m3] * m3 / m2 = [m2 / s]
  
  // Get a surface cell at (theta phi)
  Cell_t s = d.GetSurfaceCell(0.0, PIGREEK/2);
  d.PrintCell(s);


  // Get list of cells at a given Longitude
  std::vector<Cell_t> cells = d.GetCellsLongitude(PIGREEK*45.0/180.0);

  // Number of cells with non-zero activity
  int Ncells = 0;
  for(unsigned int i = 0; i < cells.size(); i++){
    if(!d.IsEqual(cells[i],s) &&
       cells[i].a238U > 0){
      Ncells++;
    }
  }
  std::cout << "Total cells to evaluate: " << Ncells << std::endl;

  // Integrate contributions from all other cells
  // Anti-Neutrino energy
  double e = 0.003; // GeV
  double FluxU238 = 0;
  double Prob = 1.0; 
  double dx, dy, dz, dr2;
  for(unsigned int i = 0; i < cells.size(); i++){
    if(!d.IsEqual(cells[i],s) &&
       cells[i].a238U > 0){
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

      std::cout << "Probability: " << Prob << std::endl;
      FluxU238 += Prob * cells[i].a238U * cells[i].rho / dr2; // [1/kg]*[kg/km3]/[km2] = [1/km5]
    }
  }
  // apply constants
  //[1/km5] * [km3] * [1] * [1/s] * = [1/(km2*s)]
  FluxU238 *= (d.m_DCell*d.m_DCell*d.m_DCell) * d.U238.Mnu * d.U238.Lamb / (4*PIGREEK);

  std::cout << "Flux U238 : " << FluxU238 << " aneutrinos / (km2*s)" << " = " << FluxU238/(100000.*100000.*1.e+06) << " aneutrinos / (cm2*us) " <<  std::endl;

  return 0;
}
