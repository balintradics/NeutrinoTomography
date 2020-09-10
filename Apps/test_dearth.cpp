#include "../EarthModel/DiscreteEarth.h"

int main(){

  DiscreteEarth d(100.0); // km cell size
  //d.PrintDensityR();
  d.SaveEarthToCSV();
  return 0;
}
