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
  DiscreteEarth d(200.0); // km cell size

  std::vector<Cell_t> scellslong = d.GetSurfaceCellsLongitude(10.0*PIGREEK/180.0);

  ofstream outfile ;
  outfile.open("Surf_Cells_Longitudet.dat");
  for(int i = 0; i < scellslong.size();i++){
    outfile << scellslong[i].x << "\t" << scellslong[i].y << "\t" << scellslong[i].z << "\t" << endl;
  }

  outfile.close();

  return 0;
}
