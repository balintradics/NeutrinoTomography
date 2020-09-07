#ifndef DISCRETEEARTH_H
#define DISCRETEEARTH_H

// goal: to have an object that can store at a given (r, phi theta)
//  - the local density
//  - the local elemental/chemical composition
//  - can plot or output the profiles of these along r, phi or theta

// build:
// $ g++ -c DiscreteEarth.cxx -I.
// $ g++ -o test_dearth test_dearth.cxx -I. DiscreteEarth.o

// ------------------------------------------------------
// Geoneutrino flux model based on O. Sramek et al,
// Earth and Planetary Science Letters 361 (2013) 356-366
// ------------------------------------------------------
//
// Antineutrino Flux model: 
// Flux(r) = (nX * LambX/4pi) * <P> * \int aX(r')*rho(r')/(|r-r'|)^2 dr'
//
// nX: number of antineutrinos per decay chain
// LambX: decay constant (1/Lifetime)
// <P> : average survival probability of neutrinos at r
// aX: Abundance of radioactive isotope (number of atoms of radioactive isotope
//     per unit mass of rock)
// rho: rock density
// [X: is a placeholder for an isotope]
//
// Note: Key issue is how to calculate the local abundance of rad. isotopes, aX
//
// ------------------------------------------------------
// aX = AX * XX / MX
//
// AX: Elemental abundance (mass of element per unit mass of rock)
// XX: Isotopic ratio (atoms of radionuclide per atoms of element)
// MX: Atomic mass
// ------------------------------------------------------
// The unknowns in the model are:
// - AX: elemental abundance
// - rho: rock density




// constants related to Earth
// from PREM model: Physics of the Earth and Planetary Interiors, 25 (1981) 297—356
const float R_E=6371.; // Earth radius in Km
const float R_LM = 3480.; // Lower Mantle / Core boundary radius in Km
const float R_Crust = 6346.6; // Crust / Mantle boundary radius in Km
const float R_Ocean = 6368.0; // Ocean / Crust boundary radius in Km




// A cell of the discretized Earth
struct Cell_t{
  float x;
  float y;
  float z;
  float rho; // density
  float a238U; // abundance of 238U rad. isotope
  float a235U; // abundance of 238U rad. isotope
  float a232Th; // abundance of 238U rad. isotope
  float a40K; // abundance of 238U rad. isotope
};
  
class DiscreteEarth {
 public:
  // constructor, parameter: dcell size of a cube in km
  DiscreteEarth(float dcell);
  ~DiscreteEarth();


  void ToSpherical(float x, float y, float z, float * r, float * theta, float * phi);
  void ToCartesian( float * theta, float * phi, float x, float y, float z);
  void SetUniformMantle();

  void SaveEarthToCSV(const char * ofilename = "earth_model.csv");
  void SaveDensityRToCSV(const char * ofilename = "density.csv");
  void SaveFluxToCSV(const char * ofilename = "surf_flux.csv");// latitude, longitude, flux

 public:
  // internally we represent the coordinates
  // with x, y, z Cartesian coordinates
  Cell_t * m_EarthCells;
  float m_DCell; // cell size
  int m_NCells; // number of cells
  
  
};


#endif // DICRETEEARTH_H
