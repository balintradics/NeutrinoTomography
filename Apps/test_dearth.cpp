#include "../EarthModel/DiscreteEarth.h"

int main(){

  DiscreteEarth d(100.0); // km cell size
  //d.PrintDensityR();
  d.SaveEarthToCSV();

  // pseudo-code
  // loop over each cell:
  //     use the properties of that cell (iosotopic activity) to decide
  //     how many neutrinos to shoot: N;
  //     shoot N  neutrinos from that cell (to a given direction);
  //     get surface surival probability at a given lat, long;

  return 0;
}
