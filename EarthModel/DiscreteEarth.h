#ifndef DISCRETEEARTH_H
#define DISCRETEEARTH_H

#include <vector>

// DiscreteEarth class: provides an object that can store at a given (r, phi theta)
//  - the local density
//  - the local elemental/chemical composition
//  - can plot or output the profiles of these along r, phi or theta
// This is to be used in conjuction with a neutrino propagation in material

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
// AX: Elemental abundance (mass of element per unit mass of rock) [dimless]
// XX: Isotopic ratio (atoms of radionuclide per atoms of element) [dimless]
// MX: Atomic mass
// ------------------------------------------------------
// The unknowns in the model are:
// - AX: elemental abundance
// - rho: rock density

// Atomic constants
struct AtomConst_t{
  float X; // Isotopic ratio
  float M; // Mass in amu = 1.661e-27 kg
  float Tau; // Half life, Gyr
  float Lamb; // Decay const, [1e-18 sec^-1]
  float Mnu; // Number of e- antineutrinos per chain
};



// constants related to Earth
// from PREM model: Physics of the Earth and Planetary Interiors, 25 (1981) 297—356
const float R_E=6371.; // Earth radius in Km
const float R_Ocean = 6368.0; // Ocean / Crust boundary radius in Km
const float R_Crust = 6346.6; // Crust / Mantle boundary radius in Km
const float R_LM = 3480.; // Lower Mantle / Core boundary radius in Km

// A cell of the discretized Earth
struct Cell_t{
  float x; // [km]
  float y; // [km]
  float z; // [km]
  float rho; // density [g/cm3]
  float a238U; // rel. abundance of 238U rad. isotope 
  float a235U; // rel. abundance of 238U rad. isotope
  float a232Th; // rel. abundance of 238U rad. isotope
  float a40K; // rel. abundance of 238U rad. isotope
};
  
// Spherical coordinate system
// r, theta (polar), phi (azimuth)
// theta: 0 - pi
// phi:  0 - 2*pi

class DiscreteEarth {
 public:
  // constructor, parameter: dcell size of a cube in km
  DiscreteEarth(float dcell);
  ~DiscreteEarth();

  void PrintCell(Cell_t cell);
  Cell_t GetRandomCell();
  Cell_t GetSurfaceCell(float theta, float phi);
  Cell_t GetCell(float x, float y, float z);
  void ToSpherical(float x, float y, float z, float * r, float * theta, float * phi);
  void ToCartesian(float r, float theta, float phi, float * x, float * y, float * z);
  void SetUniformMantle();

  void SaveEarthToCSV(const char * ofilename = "earth_model.csv");
  void SaveDensityRToCSV(const char * ofilename = "density.csv");
  void SaveFluxToCSV(const char * ofilename = "surf_flux.csv");// latitude, longitude, flux

  float SetOriginTarget(Cell_t ocell, Cell_t tcell);

  // Return the density "along" a vectorial path during propagation
  // needs to be static in order to be callable from external sources
  static double Density(double x, double y, double z);
  static double DensityAlong(double s); 

 public:
  // internally we represent the coordinates
  // with x, y, z Cartesian coordinates
  Cell_t * m_EarthCells;
  float m_DCell; // cell size
  int m_NCells; // number of cells

  // variables to aide the neutrino propagation between Earth Cells
  // Origin and Target
  Cell_t m_Tcell; // Target
  static Cell_t m_Ocell; // Origin
  static float m_Dx; // Direction x
  static float m_Dy; // Direction y
  static float m_Dz; // Direction z
  static float m_PathLength;

  // Structures to hold Atomic constants
  AtomConst_t U238;
  AtomConst_t U235;  
  AtomConst_t Th232;
  AtomConst_t K40;
  

};


#endif // DICRETEEARTH_H
