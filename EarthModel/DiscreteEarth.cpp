#include "math.h"
#include "stdlib.h"
#include <time.h> 
#include <cstdlib> 
#include <ctime>

#include "DiscreteEarth.h"

using namespace std;

DiscreteEarth::DiscreteEarth(double dcell){

  m_DCell = dcell;// [km]
  if (m_DCell < 50.0){
    cout << "Cell size is too small, to save memory resizing to 50 km" << endl;
    m_DCell = 50.0;
  }

  m_DRad = m_DCell/R_E;
  
  // How many cells we need to fill the full Earth volume?
  // We calculate it dynamically...
  m_NCells = 0;
  m_EarthCells = NULL;
  double r = 0;
  for(double ix = -R_E; ix <= R_E; ix = ix + m_DCell){
    for(double iy = -R_E; iy <= R_E; iy = iy + m_DCell){
      for(double iz = -R_E; iz <= R_E; iz = iz + m_DCell){
	// skip cells that would be outside of the Earth volume
	r = sqrt(ix*ix + iy*iy + iz*iz);
	if(r >= R_E){
	}else{
	  m_NCells++;
	}
      }
    }
  }

  
  m_EarthCells = new Cell_t[m_NCells];

  
  // Now fill the values
  // Using hard-coded density values!
  int Celli = 0;
  double rho_ret = 2.6 * 1000./1e-09; // g/cm3 -> kg/km3
  for(double ix = -R_E; ix <= R_E; ix = ix + m_DCell){
    for(double iy = -R_E; iy <= R_E; iy = iy + m_DCell){
      for(double iz = -R_E; iz <= R_E; iz = iz + m_DCell){
	r = sqrt(ix*ix + iy*iy + iz*iz);
	if(r >= R_E){
	}else{
	  
	  Cell_t cell;
	  rho_ret = Density(ix, iy,iz);
	  cell.rho = rho_ret;
	  cell.x = ix;
	  cell.y = iy;
	  cell.z = iz;
	  cell.Surf = false;
	  m_EarthCells[Celli] = cell;
	  Celli++;
	}
      }
    }
  }

  // Separately store the Surface cells where we want uniform binning in theta/phi
  m_NSurfCells = 0;
  for(double itheta = 0; itheta <= PI; itheta = itheta + m_DRad){
    for(double iphi = 0; iphi <= 2*PI; iphi = iphi + m_DRad){
      m_NSurfCells++;
    }
  }

  m_SurfCells = new Cell_t[m_NSurfCells];

  double x, y, z;
  int SurfCell = 0;
  // Separately store the Surface cells where we want uniform binning in theta/phi
  for(double itheta = 0; itheta <= PI; itheta = itheta + m_DRad){
    for(double iphi = 0; iphi <= 2*PI; iphi = iphi + m_DRad){
      Cell_t cell;
      cell.Surf = true;
      ToCartesian(R_E, itheta, iphi, &x, &y, &z);
      rho_ret = Density(x, y,z);
      cell.rho = rho_ret;
      cell.x = x;
      cell.y = y;
      cell.z = z;
      m_SurfCells[SurfCell] = cell;
      //      m_SurfCellInd.push_back(Celli);
      SurfCell++;
    }
  }


  cout << "Allocated " << m_NCells << " cells, inclusive " << m_NSurfCells << " surface cells" << endl;

  // Set Atomic constants
  U238.X = 0.9927;
  U238.M = 238.051*amu;
  U238.Tau = 4.468;
  U238.Lamb = 4.916*lu;
  U238.Mnu = 6;
  
  U235.X = 0.007204;
  U235.M =235.044*amu;
  U235.Tau = 0.704;
  U235.Lamb = 31.2*lu;
  U235.Mnu = 4;
  
  Th232.X = 1.0;
  Th232.M = 232.038*amu;
  Th232.Tau = 14.05;
  Th232.Lamb = 1.563*lu;
  Th232.Mnu = 4;
  
  K40.X = 117.0e-06;
  K40.M = 39.9640*amu;
  K40.Tau = 1.265;
  K40.Lamb = 17.36*lu;
  K40.Mnu = 0.8928;

  //Depleted Mantle composition from Arevalo and McDonough
  DepMantle.A_U = 8*ppb;// +/- 2 ppb
  DepMantle.A_Th = 22*ppb;// +/- 4 ppb
  DepMantle.A_K = 152*ppm;// +/- 30 ppm
  
  // set initial seed value to system clock
  srand( time(NULL) );
  
  // Initialize number of detector cells and cell size
  m_DetCellSize = m_DCell;
  m_NDetCells = (int)(2*R_E/m_DetCellSize);



}

DiscreteEarth::~DiscreteEarth(){

  delete[] m_EarthCells; 

}

// Get the normal vector to the plane of the Longitude
// and return it conveniently packged as a Quaternion
Quat4d_t DiscreteEarth::GetNormalToLongitude(double phi){
  // Let vec(P0) = (P0x, P0y, P0z) be a point given in the plane of the longitude
  // let vec(n) = (nx, ny, nz) an orthogonal vector to this plane
  // then vec(P) = (Px, Py, Pz) will be in the plane if (vec(P) - vec(P0)) * vec(n) = 0
  
  // We pick 2 vectors in the plane
  double P0x, P0y, P0z; // given by the longitude
  ToCartesian(R_E, 2.0, phi, &P0x,&P0y, &P0z);
  double P1x, P1y, P1z;
  ToCartesian(R_E, 2.5, phi, &P1x,&P1y, &P1z);
  double P2x, P2y, P2z;
  ToCartesian(R_E, 3.0, phi, &P2x,&P2y, &P2z);

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

  Quat4d_t normQ;
  normQ.x = nx;
  normQ.y = ny;
  normQ.z = nz;
  normQ.w = 0.0;

  return normQ;

}


void DiscreteEarth::GetLongitudePlaneBasis(double phi, Quat4d_t * basis1, Quat4d_t * basis2){

  // Rotation by Quaternion algebra:
  // rotate a vector "p" by angle "angle" around a (unit) axis "r"
  // Form Quaterniin: q1 = (px, py, pz, 0)
  // Form unit rotation Quaternion: q2 = (rx*sin(angle/2), ry*sin(angle/2), rz*sin(angle/2), cos(angle/2))
  // Forward rotated Quaternion: Q3 = q2 * q1 * q2^{*}  
  // Backward rotated Quaternion: Q3 = q2^{*} * q1 * q2 

  // Rotation axis is z - axis if rotate only by Longitude
  // Rotation Quaternion:
  Quat4d_t qr = ToRotQuaternion(0, 0, 1, phi);
  //  PrintQ(qr);

  // conjugate
  Quat4d_t qrnc = ConjugateQ(qr);
  //  PrintQ(qrnc);

  // We want to rotate the x axes into the longitude plane
  Quat4d_t qe = ToQuaternion(1.0, 0.0, 0.0);

  // do the rotation
  Quat4d_t q = MultiplyQ(qr, qe);
  Quat4d_t q3 = MultiplyQ(q, qrnc);

  // The rotated x - basis vector
  (*basis1).x = q3.x;
  (*basis1).y = q3.y;
  (*basis1).z = q3.z;
  (*basis1).w = 0.0;

  // This should be just the z - basis vector
  // The z - axis stays unchanged when rotating around it by phi
  (*basis2).x = 0;
  (*basis2).y = 0;
  (*basis2).z = 1;
  (*basis2).w = 0.0;

  return;

}

void DiscreteEarth::CreateDetector(Quat4d_t basis1, Quat4d_t basis2, Quat4d_t qrot, double angle){
  // Basis1 is the negative normal of the detector plane
  // Equation of a line: a point on the line and a vector along the line
  // Basis1*R_E gives the center point of the line
  // Construct a vector along the line by rotating the Basis1 by 90deg

  m_Det1.clear();
  m_Det2.clear();

  Quat4d_t linevec = RotateQuaternion(basis1, qrot, 90.0*PI/180.0);

  // Add half the detector cells + and - from the central point

  // note we start with 1 because we already added the center det cell
  for(int i = m_NDetCells/2 + 1; i >= 1; i--){
     Quat4d_t q;
     q.x = basis1.x*R_E - i*m_DetCellSize*linevec.x;
     q.y = basis1.y*R_E - i*m_DetCellSize*linevec.y;
     q.z = basis1.z*R_E - i*m_DetCellSize*linevec.z;
     q.w = 0.0;

     //     RotateQuaternion(&q, qrot, angle);
     
     m_Det1.push_back(q);
  }

  
  for(int i = 0; i <= m_NDetCells/2 + 1; i++){
     Quat4d_t q;
     q.x = basis1.x*R_E + i*m_DetCellSize*linevec.x;
     q.y = basis1.y*R_E + i*m_DetCellSize*linevec.y;
     q.z = basis1.z*R_E + i*m_DetCellSize*linevec.z;
     q.w = 0.0;

     //     RotateQuaternion(&q, qrot, angle);
     m_Det1.push_back(q);
  }

  // Now add fake detectors at the opposite side
  // note we start with 1 because we already added the center det cell
  for(int i = m_NDetCells/2 + 1; i >= 1; i--){
     Quat4d_t q;
     q.x = -1.0*basis1.x*R_E - i*m_DetCellSize*linevec.x;
     q.y = -1.0*basis1.y*R_E - i*m_DetCellSize*linevec.y;
     q.z = -1.0*basis1.z*R_E - i*m_DetCellSize*linevec.z;
     q.w = 0.0;

     //     RotateQuaternion(&q, qrot, angle);
     
     m_Det2.push_back(q);
  }

  
  for(int i = 0; i <= m_NDetCells/2 + 1; i++){
     Quat4d_t q;
     q.x = -1.0*basis1.x*R_E + i*m_DetCellSize*linevec.x;
     q.y = -1.0*basis1.y*R_E + i*m_DetCellSize*linevec.y;
     q.z = -1.0*basis1.z*R_E + i*m_DetCellSize*linevec.z;
     q.w = 0.0;

     //     RotateQuaternion(&q, qrot, angle);
     m_Det2.push_back(q);
  }
  


}


void DiscreteEarth::PrintDetector(){

  cout << "Number of detector cells: " << m_NDetCells << endl;
  cout << "Detector cell size: " << m_DetCellSize << endl;
  
  for(int i = 0; i < m_Det1.size(); i++){
    PrintQ(m_Det1[i]);
  }

}

// Only for Det1 - Det2 is the mirror "fake" detector
Quat4d_t DiscreteEarth::GetDetCoord(Cell_t surfcell, Quat4d_t basis2){

  double z_coord = surfcell.z * basis2.z;
  Quat4d_t q;
  bool found = false;
  // find the corresponding detector bin coordinate
  for(int i = 0; i < m_Det1.size()-1; i++){
    if(z_coord >= m_Det1[i].z && z_coord <= m_Det1[i+1].z){
      q = m_Det1[i];
      found = true;
      break;
    }
  }
  
  if(!found){
    cout << "Error: No corresponding detector coordinate!" << surfcell.x << ", " << surfcell.y << ", " << surfcell.z << endl;
  }

  return q;

}

// Get cos of angular overlap between to 3D vectors
// for parallel vectros it gives cos(theta) =  1.0
double DiscreteEarth::GetVecOverlap(Quat4d_t v1, Quat4d_t v2){

  // Collimator approximation: only accept neutrinos that are parallel
  // with the detector plane basis
  // u*v = |u|*|v|*cos(theta)
  double prod = v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
  double mag1 = sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);
  double mag2 = sqrt(v2.x*v2.x + v2.y*v2.y + v2.z*v2.z);
  double cos_theta = prod/(mag1*mag2);
  
  return cos_theta;

}


// Only for Det1 - Det2 is the mirror "fake" detector
int DiscreteEarth::GetDetBin(Cell_t surfcell, Quat4d_t basis2){

  // Project the surfcell onto the basis2 vector to get the coordinate
  // Then find the cell with that coordinate
  //  cout << "GetDetBin: Surface cell: " << surfcell.x << ", " << surfcell.y << ", " << surfcell.z << endl;
  double surf_projcoord = surfcell.x*basis2.x + surfcell.y*basis2.y + surfcell.z*basis2.z;
  
  int ret = 0;
  bool found = false;
  double diff = 99999999.0;
  // find the corresponding detector bin coordinate
  for(int i = 0; i < m_Det1.size()-1; i++){
    double det_projcoord = m_Det1[i].x*basis2.x + m_Det1[i].y*basis2.y + m_Det1[i].z*basis2.z;
    
    diff = fabs(det_projcoord - surf_projcoord);
    //    cout << "GetDetBin: Det[" << i << "], current diff: " << surf_projcoord - det_projcoord << endl;
    if(diff  < m_DCell){
      ret = i;
      found = true;
      //      cout << "Found: diff = " << diff << endl;
      break;
    }
  }
  
  if(!found){
    cout << "Error: No corresponding detector coordinate!" << surfcell.x << ", " << surfcell.y << ", " << surfcell.z << ", diff: " << diff << endl;
    
  }

  return ret;

}




// Outputs the list of cells lying in the plane of a longitude
void DiscreteEarth::SaveCellsLongitudeToFile(double phi_fix, const char * ofilename){
  
  Quat4d_t normQ = GetNormalToLongitude(phi_fix);

  // A vector lying in the plane of the longitude
  double P0x, P0y, P0z; // given by the longitude radial vector
  ToCartesian(R_E, 2.0, phi_fix, &P0x,&P0y, &P0z);

  double dx, dy, dz, dMag;
  double prod = 0;
  ofstream outfile;
  outfile.open(ofilename);
  for(unsigned int i = 0; i < m_NCells; i++){
    // only save cells that lie in the (vicinity of the) plane
    dx = (m_EarthCells[i].x - P0x);
    dy = (m_EarthCells[i].y - P0y);
    dz = (m_EarthCells[i].z - P0z);
    dMag = sqrt(normQ.x*normQ.x + normQ.y*normQ.y + normQ.z*normQ.z);
    dx /= dMag; dy /= dMag; dz /= dMag;
    prod = normQ.x*dx + normQ.y*dy + normQ.z*dz; // inner product between normal veector of plane and current cell
    //    std::cout << m_EarthCells[i].x << ", " << m_EarthCells[i].y << ", " << m_EarthCells[i].z << ": " << prod << std::endl;
    if( fabs(prod) <= 100 ){
      outfile  << m_EarthCells[i].x << "\t" << m_EarthCells[i].y << "\t" << m_EarthCells[i].z << "\t" << m_EarthCells[i].a238U <<  endl;
    }
  }
  outfile.close();

}


//Outputs the list of cells in the Y-Z plane
void DiscreteEarth::SaveActivityMap2DToFile(const char * ofilename){
    
  ofstream outfile;
  outfile.open(ofilename);
  double r = 0;
  for(unsigned int i = 0; i < m_NCells; i++){
    // only save along x = 0
    if( fabs(m_EarthCells[i].x) <= m_DCell){
      outfile  << m_EarthCells[i].y << "\t" << m_EarthCells[i].z << "\t" << m_EarthCells[i].a238U <<  endl;
    }
  }
  outfile.close();

}

//Outputs the list of cells in the Y-Z plane
void DiscreteEarth::SaveSurfaceCellsToFile(const char * ofilename){
    
  ofstream outfile;
  outfile.open(ofilename);
  for(unsigned int i = 0; i < m_NSurfCells; i++){
    outfile  << m_SurfCells[i].x << "\t" << m_SurfCells[i].y << "\t" << m_SurfCells[i].z << "\t" << m_SurfCells[i].a238U <<  endl;
  }
  outfile.close();

}




void DiscreteEarth::SaveDensityRToFile(const char * ofilename){
    
  ofstream outfile;
  outfile.open(ofilename);
  double r = 0;
  for(unsigned int i = 0; i < m_NCells; i++){
    // only print along x = 0, y = 0
    if( fabs(m_EarthCells[i].x) <= m_DCell && fabs(m_EarthCells[i].y) <= m_DCell && m_EarthCells[i].z >= 0){
      r = sqrt(m_EarthCells[i].x*m_EarthCells[i].x+m_EarthCells[i].y*m_EarthCells[i].y+m_EarthCells[i].z*m_EarthCells[i].z);
      outfile  << m_EarthCells[i].z << "\t" << m_EarthCells[i].rho << endl;
    }
  }
  outfile.close();

}

void DiscreteEarth::SaveEarthToFile(const char * ofilename){
  cout << "Saving : " << m_NCells << " cells values to " << ofilename << endl;    
  ofstream outfile;
  outfile.open(ofilename);
  for(unsigned int i = 0; i < m_NCells; i++){
    outfile  << m_EarthCells[i].x << ", " << m_EarthCells[i].y << ", " << m_EarthCells[i].z << ", " << m_EarthCells[i].rho << endl;
  }
  outfile.close();

}


// Example how to save output file in format: longitude  latitude , the flux/other value
void DiscreteEarth::SaveFluxMap(const char * ofluxfilename){
  cout << "Saving values to fixed coordinate (lat/lon) map: " << ofluxfilename << endl;

  // Latitude (theta): -pi/2 - pi/2 (-90 - 90 degrees)
  // Longitude (phi): 0 - 2*pi (0 - 360 degrees)

  // Note: correspondence between theta and latitude: latitude = theta - pi/2

  ofstream outfile;
  outfile.open(ofluxfilename);
  double x, y, z, r, theta, phi;

  for(int i = 0; i < m_NSurfCells;i++){
    Cell_t cell = m_SurfCells[i];
    ToSpherical(cell.x, cell.y, cell.z, &r, &theta, &phi );
    outfile << -theta + PI/2 << "\t" << phi << "\t" << cell.a238U << endl;
  }

  outfile.close();

}

void DiscreteEarth::ToSpherical(double x, double y, double z, double * r, double * theta, double * phi){
  *r = sqrt(x*x + y*y + z*z);
  *theta = acos(z/(*r));// returns [0, pi]
  *phi = atan2(y,x) + PI; // returns [0,2pi]; atan2 alone returns [-pi, pi]
}

void DiscreteEarth::ToCartesian(double r, double theta, double phi, double * x, double * y, double * z){
  *x = r*sin(theta)*cos(phi);
  *y = r*sin(theta)*sin(phi);
  *z = r*cos(theta);
}



void DiscreteEarth::PrintCell(Cell_t cell){
  cout << "-----------" << endl;
  cout << " E A R T H  C E L L: " << endl;
  cout << "-----------" << endl;
  cout << " X= " <<  cell.x << " km " << endl;
  cout << " Y= " <<  cell.y << " km " << endl;
  cout << " Z= " <<  cell.z << " km " << endl;
  if(cell.Surf){
    cout << " Surface cell = " << " True " << std::endl;
  }else{
    cout << " Surface cell = " << " False " << std::endl;
  }
  cout << " rho= " << cell.rho << " kg/km3 " << endl;
  cout << " a_238U= " << cell.a238U << " /kg " << endl;
  cout << " a_235U= " << cell.a235U << " /kg " << endl;
  cout << " a_232Th= " << cell.a232Th << " /kg" << endl;
  cout << " a_40K= " << cell.a40K << " /kg " << endl;

  cout << "-----------" << endl;
}

Cell_t DiscreteEarth::GetRandomCell(){
  /* initialize random seed: */
  srand (time(NULL));
  return m_EarthCells[rand() % m_NCells];

}

Cell_t DiscreteEarth::GetSurfaceCell(double theta, double phi){
  double x, y, z;
  ToCartesian(R_E, theta, phi, &x, &y, &z);

  //cout << "Looking for surface cell at: (r, theta, phi): " << R_E << ", " << theta << ", " << phi << endl;
  //cout << "(x, y, z): " << x << ", " << y << ", " << z << endl;
  
  Cell_t result;
  bool found = false;
  for(unsigned int i = 0; i < m_NSurfCells; i++){
    if( fabs(m_SurfCells[i].x-x) <= m_DCell && fabs(m_SurfCells[i].y-y) <= m_DCell && fabs(m_SurfCells[i].z-z) <= m_DCell){
      result = m_SurfCells[i];
      found = true;
      break;
    }
  }
  if(!found){
    cout<< "Cell not found" << endl;
  }
  return result;
}

Cell_t DiscreteEarth::GetCell(double x, double y, double z){
  //cout << "Looking for surface cell at: (x, y, z): " << x << ", " << y << ", " << z << endl;
  
  Cell_t result;
  bool found = false;
  for(unsigned int i = 0; i < m_NCells; i++){
    // only print along x = 0, y = 0
    if( fabs(m_EarthCells[i].x-x) <= m_DCell && fabs(m_EarthCells[i].y-y) <= m_DCell && fabs(m_EarthCells[i].z-z) <= m_DCell){
      result = m_EarthCells[i];
      found = true;
      break;
    }
  }
  if(!found){
    cout<< "Cell not found" << endl;
  }
  return result;

}


// Returns the list of cells lying in the plane of a longitude
std::vector<Cell_t> DiscreteEarth::GetCellsLongitude(double phi_fix){
  // Let vec(P0) = (P0x, P0y, P0z) be a point given in the plane of the longitude
  // let vec(n) = (nx, ny, nz) an orthogonal vector to this plane
  // then vec(P) = (Px, Py, Pz) will be in the plane if (vec(P) - vec(P0)) * vec(n) = 0
  
  // We pick 2 vectors in the plane
  double P0x, P0y, P0z; // given by the longitude
  ToCartesian(R_E, 2.0, phi_fix, &P0x,&P0y, &P0z);
  double P1x, P1y, P1z;
  ToCartesian(R_E, 2.5, phi_fix, &P1x,&P1y, &P1z);
  double P2x, P2y, P2z;
  ToCartesian(R_E, 3.0, phi_fix, &P2x,&P2y, &P2z);

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
  
  double dx, dy, dz, dMag;
  double prod = 0;

  std::vector<Cell_t> cell_vec;
  for(unsigned int i = 0; i < m_NCells; i++){
    // only save cells that lie in the (vicinity of the) plane
    dx = (m_EarthCells[i].x - P0x);
    dy = (m_EarthCells[i].y - P0y);
    dz = (m_EarthCells[i].z - P0z);
    dMag = sqrt(nx*nx + ny*ny + nz*nz);
    dx /= dMag; dy /= dMag; dz /= dMag;
    prod = nx*dx + ny*dy + nz*dz;

    if( fabs(prod) <= 100 ){
      cell_vec.push_back(m_EarthCells[i]);
    }
  }
  
  return cell_vec;

}


std::vector<Cell_t> DiscreteEarth::GetSurfaceCellsLongitude(double phi_fix){
  // Let vec(P0) = (P0x, P0y, P0z) be a point given in the plane of the longitude
  // let vec(n) = (nx, ny, nz) an orthogonal vector to this plane
  // then vec(P) = (Px, Py, Pz) will be in the plane if (vec(P) - vec(P0)) * vec(n) = 0
  
  std::vector<Cell_t> cell_vec;
  double r, theta, phi;
  for(unsigned int i = 0; i < m_NSurfCells; i++){
    // only save cells that lie in the (vicinity of the) plane
    ToSpherical(m_SurfCells[i].x, m_SurfCells[i].y, m_SurfCells[i].z, &r, &theta, &phi);
    if( fabs(phi - phi_fix) < m_DRad/2.0 || fabs(phi - (phi_fix+PI)) < m_DRad/2.0){
      cell_vec.push_back(m_SurfCells[i]);
    }
  }
  
  return cell_vec;



}


bool DiscreteEarth::IsEqual(Cell_t a, Cell_t b) const{
  if(a.x == b.x && a.y == b.y && a.z == b.z){
    return true;
  }
  return false;

}

double DiscreteEarth::SetOriginTarget(Cell_t ocell, Cell_t tcell){
  m_Ocell = ocell;
  m_Tcell = tcell;

 // Calculate the difference vector pointing
 // to the surface cell from the original cell
  m_Dx = m_Tcell.x - m_Ocell.x;
  m_Dy = m_Tcell.y - m_Ocell.y;
  m_Dz = m_Tcell.z - m_Ocell.z;

  // Set the length
  m_PathLength = sqrt(m_Dx*m_Dx + m_Dy*m_Dy + m_Dz*m_Dz);

  return m_PathLength;

}


double DiscreteEarth::Density(double x, double y, double z){
  // Calculate fast the density at this point
  double rho_ret = 2.6 * 1000./1.e-09;// g/cm3 -> kg/km3

  double r = sqrt(x*x+y*y+z*z);
  double R =0;
  if(r > R_E){
    rho_ret = 0.0;
  }else{
    R = r/R_E;
    if(r>0. && r<1221.5) 
      rho_ret = (13.0885-8.8381*R*R);
    else if (r >=1221.5 && r<3480.)
      rho_ret = (12.5815-1.2638*R-3.6426*R*R-5.5281*R*R*R);
    else if (r>=3480. && r<5701.)
      rho_ret = (7.9565-6.4761*R + 5.5283*R*R-3.0807*R*R*R);
    else if (r >=5701.0 && r<5771.0)
      rho_ret = (5.3197-1.4836*R);
    else if (r>=5771.0 && r<5971.0)
      rho_ret = (11.2494-8.0298*R);
    else if (r>=5971.0 && r<6151.0)
      rho_ret = (7.1089-3.8045*R);
    else if (r>=6151.0 && r<6346.6)
      rho_ret = (2.6910+0.6924*R);
    else if (r>=6346.6 && r<6356.6)
      rho_ret = 2.9;
    else if (r>=6356.6 && r <=R_E)
      rho_ret = 2.6;
  }

  return rho_ret * 1000./1.e-09; // kg/km3

}

// Return the density "along" a vectorial path during propagation
// It is used in the neutrino oscillation integrator on the r.h.s
// of the diff. equation
double DiscreteEarth::DensityAlong(double s){
  // The input "s" is the position along the path 
  // We just have to find the corresponding Cell and return the local
  // density inside

  double fractionalPath = s/m_PathLength;
  // Current location is vOrigin + fractionalPath * vDirection
  double m_Px = m_Ocell.x+fractionalPath*m_Dx;
  double m_Py = m_Ocell.y+fractionalPath*m_Dy;
  double m_Pz = m_Ocell.z+fractionalPath*m_Dz;

  // Note: NeutrinoOscillation part needs densities in g/cm3!
  double rho_ret = Density(m_Px, m_Py, m_Pz)/( 1000.*1e+09); // kg/km3 -> g/cm3

  //  cout << fractionalPath << ": " <<  m_Px << "\t" << m_Py << "\t" << m_Pz << "\t" << rho_ret << endl;

  return rho_ret;

}

void DiscreteEarth::SetUniformMantle(Comp_t comp){

  // Simplest model sets the Cell radiogenic compositions
  // to uniform in the mantle

  double r = 0;
  for(unsigned int i = 0; i < m_NCells; i++){
    r = sqrt(m_EarthCells[i].x*m_EarthCells[i].x + m_EarthCells[i].y*m_EarthCells[i].y + m_EarthCells[i].z*m_EarthCells[i].z);
    if(r > R_LM){
      m_EarthCells[i].a238U = comp.A_U*U238.X/U238.M;
      m_EarthCells[i].a235U = comp.A_U*U235.X/U235.M;
      m_EarthCells[i].a232Th = comp.A_Th*Th232.X/Th232.M;
      m_EarthCells[i].a40K = comp.A_K*K40.X/K40.M;
      //      cout << "R: " << r << "\t " << m_EarthCells[i].a238U << endl;
    }
  }

}


// An idealized, single 1000 km thick pile 
// with a lateral extent: Lat = 0 - 76 degrees
void DiscreteEarth::SetMantleP1(){
  // Enrichment factor E_X = A_X^{EM} / A_X^{DM}
  // where A_X^{DM} is the normal Depleted Mantle composition
  // A_X^{EM} is the enriched Mantle composition
  
  //Enriched Mantel composition from Arevalo and McDonough

  double E_X_U = 6.3;
  double E_X_Th = 12;

  // Add first a uniform mantle 
  SetUniformMantle(DepMantle);

  // Overwrite the uniform mantle where enrichment is needed
  Comp_t EnMantle = DepMantle;
  EnMantle.A_U *= E_X_U;
  EnMantle.A_Th *= E_X_Th;
  SetMantleEnrichedLayer(EnMantle, 1000.0, 90.0*PI/180, 120*PI/180, 0.0*PI/180.0, 20*PI/180.0);
  // SetMantleEnrichedLayer(EnMantle, 1000.0, 90.0*PI/180, 120*PI/180, 0.0*PI/180.0, 20*PI/180.0);

  // In the internal coord. system: theta = 0 corresponds to upwards (in geography that's +90 Latitude)
}



// An idealized, double 1000 km thick piles 
// with a lateral extent: Lat = 0 - 76 degrees at both sides
void DiscreteEarth::SetMantleP2(){
  // Enrichment factor E_X = A_X^{EM} / A_X^{DM}
  // where A_X^{DM} is the normal Depleted Mantle composition
  // A_X^{EM} is the enriched Mantle composition
  
  //Enriched Mantel composition from Arevalo and McDonough

  double E_X_U = 6.3;
  double E_X_Th = 12;

  // Add first a uniform mantle 
  SetUniformMantle(DepMantle);

  // Overwrite the uniform mantle where enrichment is needed
  Comp_t EnMantle = DepMantle;
  EnMantle.A_U *= E_X_U;
  EnMantle.A_Th *= E_X_Th;
  SetMantleEnrichedLayer(EnMantle, 1000.0, 60.0*PI/180, 120*PI/180, 0.0*PI/180.0, 20*PI/180.0);
  SetMantleEnrichedLayer(EnMantle, 1000.0, 60.0*PI/180, 120*PI/180, (0.0+180.0)*PI/180.0, (20+180)*PI/180.0);

  // To set North/South enriched layers:
  // SetMantleEnrichedLayer(EnMantle, 1000.0, (0.0)*PI/180, (30)*PI/180, 0.0*PI/180.0, 20*PI/180.0);
  // SetMantleEnrichedLayer(EnMantle, 1000.0, (0.0)*PI/180, (30)*PI/180, (0.0+180)*PI/180.0, (20+180)*PI/180.0);

  // SetMantleEnrichedLayer(EnMantle, 1000.0, (150.0)*PI/180, (180.0)*PI/180, 0.0*PI/180.0, 20*PI/180.0);
  // SetMantleEnrichedLayer(EnMantle, 1000.0, (150.0)*PI/180, (180.0)*PI/180, (0.0+180)*PI/180.0, (20+180)*PI/180.0);

  // In the internal coord. system: theta = 0 corresponds to upwards (in geography that's +90 Latitude)
}


// Set single enriched layer
void DiscreteEarth::SetMantleEnrichedLayer(Comp_t comp, double deltaR, double latMin, double latMax, double lonMin, double lonMax ){

  double r = 0;
  double phi, theta;
  for(unsigned int i = 0; i < m_NCells; i++){
    ToSpherical(m_EarthCells[i].x, m_EarthCells[i].y, m_EarthCells[i].z, &r, &theta, &phi );
    if(r >= R_LM && r <= (R_LM + deltaR) && 
       (theta  )>= latMin && (theta ) <= latMax &&
       phi >= lonMin && phi <= lonMax ){
      m_EarthCells[i].a238U = comp.A_U*U238.X/U238.M;
      m_EarthCells[i].a235U = comp.A_U*U235.X/U235.M;
      m_EarthCells[i].a232Th = comp.A_Th*Th232.X/Th232.M;
      m_EarthCells[i].a40K = comp.A_K*K40.X/K40.M;
      //      cout << "R, Lat, long: " << r << "\t " << theta - PI/2 << "\t" << phi << "\t" << m_EarthCells[i].a238U << endl;

    }

  }

}

// Quaternion algebra for 3D rotation about arbitrary axis
Quat4d_t DiscreteEarth::ToQuaternion(double x, double y, double z){
  Quat4d_t q;
  q.x = x;
  q.y = y;
  q.z = z;
  q.w = 0.0;
  return q;
}


// Rotation quaternion should be based on normalized components
Quat4d_t DiscreteEarth::ToRotQuaternion(double x, double y, double z, double angle){
  Quat4d_t q;
  double s = sin(angle/2.0);
  q.x = x * s;
  q.y = y * s;
  q.z = z * s;
  q.w = cos(angle/2);
  return q;
}


Quat4d_t DiscreteEarth::ConjugateQ(Quat4d_t q1){
  Quat4d_t qc;
  qc.x = -q1.x;
  qc.y = -q1.y;
  qc.z = -q1.z;
  qc.w = q1.w;
  return qc;

}

Quat4d_t DiscreteEarth::NormaliseQ(Quat4d_t q1){
  Quat4d_t qn;
  double n = sqrt(q1.x*q1.x + q1.y*q1.y + q1.z*q1.z);
  qn.x = q1.x / n;
  qn.y = q1.y / n;
  qn.z = q1.z / n;
  qn.w = q1.w / n;
  return qn;
}


Quat4d_t DiscreteEarth::ScaleQ(Quat4d_t q1, double s){
  Quat4d_t qs;
  qs.x = q1.x * s;
  qs.y = q1.y * s;
  qs.z = q1.z * s;
  qs.w = q1.w * s;
  return qs;

}

Quat4d_t DiscreteEarth::MultiplyQ(Quat4d_t q1, Quat4d_t q2){
  Quat4d_t qm;
  qm.x =  q1.x * q2.w + q1.y * q2.z - q1.z * q2.y + q1.w * q2.x;
  qm.y = -q1.x * q2.z + q1.y * q2.w + q1.z * q2.x + q1.w * q2.y;
  qm.z =  q1.x * q2.y - q1.y * q2.x + q1.z * q2.w + q1.w * q2.z;
  qm.w = -q1.x * q2.x - q1.y * q2.y - q1.z * q2.z + q1.w * q2.w;
  return qm;
}

void DiscreteEarth::PrintQ(Quat4d_t q1){
  std::cout << "Quaternion: " << std::endl;
  std::cout << "\t x: " << q1.x << std::endl;
  std::cout << "\t y: " << q1.y << std::endl;
  std::cout << "\t z: " << q1.z << std::endl;
  std::cout << "\t w: " << q1.w << std::endl;
  std::cout << "\t Norm: " << sqrt(q1.x*q1.x+q1.y*q1.y+q1.z*q1.z + q1.w*q1.w) << std::endl;
}

// Rotate a single quaternion q around the rotation axis given by rotq by angle
Quat4d_t DiscreteEarth::RotateQuaternion(Quat4d_t quat, Quat4d_t rotq, double angle){

  // Rotation by Quaternion algebra:
  // rotate a vector "p" by angle "angle" around a (unit) axis "r"
  // Form Quaterniin: q1 = (px, py, pz, 0)
  // Form unit rotation Quaternion: q2 = (rx*sin(angle/2), ry*sin(angle/2), rz*sin(angle/2), cos(angle/2))
  // Forward rotated Quaternion: Q3 = q2 * q1 * q2^{*}  
  // Backward rotated Quaternion: Q3 = q2^{*} * q1 * q2 

  // Rotation Quaternion:
  Quat4d_t qr = ToRotQuaternion(rotq.x, rotq.y, rotq.z, angle);
  // PrintQ(qr);

  // conjugate
  Quat4d_t qrnc = ConjugateQ(qr);
  //  PrintQ(qrnc);


  // Form Quaternion
  Quat4d_t qe = quat;
  // Rotate
  Quat4d_t q = MultiplyQ(qr, qe);
  Quat4d_t qfinal = MultiplyQ(q, qrnc);

  return qfinal;
  

}

  



// Rotate all cells of Earth by theta angle around an axis, rotation vector (rx, ry, rz) must be normalized
void DiscreteEarth::RotateEarth(double angle, double rx, double ry, double rz){

  // Rotation by Quaternion algebra:
  // rotate a vector "p" by angle "angle" around a (unit) axis "r"
  // Form Quaterniin: q1 = (px, py, pz, 0)
  // Form unit rotation Quaternion: q2 = (rx*sin(angle/2), ry*sin(angle/2), rz*sin(angle/2), cos(angle/2))
  // Forward rotated Quaternion: Q3 = q2 * q1 * q2^{*}  
  // Backward rotated Quaternion: Q3 = q2^{*} * q1 * q2 

  // Rotation Quaternion:
  Quat4d_t qr = ToRotQuaternion(rx, ry, rz, angle);
  // PrintQ(qr);

  // conjugate
  Quat4d_t qrnc = ConjugateQ(qr);
  //  PrintQ(qrnc);

  // Loop through all Earth Cells 
  for(int i = 0; i < m_NCells; i++){
    // Form Quaternion
    Quat4d_t qe = ToQuaternion(m_EarthCells[i].x, m_EarthCells[i].y, m_EarthCells[i].z);
    // Rotate
    Quat4d_t q = MultiplyQ(qr, qe);
    Quat4d_t q3 = MultiplyQ(q, qrnc);
    // Set the new coordinates of the Earth cells
    m_EarthCells[i].x = q3.x;
    m_EarthCells[i].y = q3.y;
    m_EarthCells[i].z = q3.z;
  }

  // Loop through all Surface Cells
  for(int i = 0; i < m_NSurfCells; i++){
    // Form Quaternion
    Quat4d_t qe = ToQuaternion(m_SurfCells[i].x, m_SurfCells[i].y, m_SurfCells[i].z);
    // Rotate
    Quat4d_t q = MultiplyQ(qr, qe);
    Quat4d_t q3 = MultiplyQ(q, qrnc);
    // Set the new coordinates of the Surf cells
    m_SurfCells[i].x = q3.x;
    m_SurfCells[i].y = q3.y;
    m_SurfCells[i].z = q3.z;
  }


}

void DiscreteEarth::ReadOscProb(const char * infilename){

  ifstream infile;
  infile.open(infilename);
  double x, y, z, theta, phi, prob;
  while(1){
    infile >> x >> y >> z >> theta >> phi >> prob;
    if(!infile.good())break;
    m_Px.push_back(x);
    m_Py.push_back(y);
    m_Pz.push_back(z);
    m_Ptheta.push_back(theta);
    m_Pphi.push_back(phi);
    m_Prob.push_back(prob);
  }

  cout << "Read " << m_Px.size() << " probabilty values from file "  << infilename << endl;


}




// Get oscillation probability between surface cell at {theta, phi} and an Earth Cell {x, y, z}
// if it is not found then do calculate it
double DiscreteEarth::GetOscProb(double theta_surf, double phi_surf, double x_cell, double y_cell, double z_cell){
  
  double osc_prob = 0.544; // mean probability from Earth and Planetary Science Letters 361 (2013) 356-366
  double eps = m_DCell/2.0;
  double epsr = m_DRad/2.0;
  bool found = false;
  // Find nearest neighbour...
  double x, y, z, theta, phi;
  for(unsigned int i = 0; i < m_Px.size(); i++){
    
    if ( fabs(m_Px[i] - x_cell) < eps && 
	 fabs(m_Py[i] - y_cell) < eps && 
	 fabs(m_Pz[i] - z_cell) < eps && 
	 fabs(m_Ptheta[i] - theta_surf) < epsr  && 
	 fabs(m_Pphi[i] - phi_surf) < epsr ){
      osc_prob = m_Prob[i];
      x = m_Px[i];
      y = m_Py[i];
      z = m_Pz[i];
      theta = m_Ptheta[i];
      phi = m_Pphi[i];
      found = true;
      break;
    }
  }

  if(found){
    cout << "Yaayy! : "  << x << ", " << y << ", " << z << ";;; " << theta << ", " << phi << " ;;; " << osc_prob << endl;
    return osc_prob;
  }else{
    return -1;
  }

  return osc_prob;


}

void DiscreteEarth::WriteOscProb(double theta_surf, double phi_surf, double x_cell, double y_cell, double z_cell, double prob){

  cout << "Writing source/target/prob: " << x_cell << ", " << y_cell << ", " << z_cell << ", " << theta_surf << ", " << phi_surf << " ;;; " << prob << endl;
  m_Px.push_back(x_cell);
  m_Py.push_back(y_cell);
  m_Pz.push_back(z_cell);
  m_Ptheta.push_back(theta_surf);
  m_Pphi.push_back(phi_surf);
  m_Prob.push_back(prob);
  
}

double DiscreteEarth::GetRandomNumber(){

  double r = ((double) rand() / (RAND_MAX));

  return r;

}

void DiscreteEarth::OpenLMF(const std::string filename){

  m_LMFile = std::make_shared<std::ofstream>(filename, std::ios::out | std::ios::binary);

  if(!(m_LMFile->is_open())) {
    std::cerr << "DiscreteEarth::OpenLMF; cannot open output file: " << filename << std::endl;
    return;
  }

  m_LMFile->clear(); //clear error bits, just in case
  const std::string header0 = "SAFIR_TOF_coincListMode";
  m_LMFile->write(&header0[0], header0.size());
  
  const std::string fill(64, 0);
  m_LMFile->write(&fill[0], fill.size());

}

void DiscreteEarth::ProcessLMF(int crys0, int crys1){

  // Process coincidence

  // const int64_t coincTimePS = coinc->getTime();
  // std::vector<std::shared_ptr<CorrectedHit> > hits = coinc->getHits();
  // if(hits.size() != 2)
  //   return;
  
  // if((coincTimePS - m_lastTime) > 1e9) { //-> one new timestamp about every ms
  //   m_lastTime = coincTimePS;
  //   ListModeData dataSetTime;
  //   dataSetTime.timeEntry.type = 1;
  //   dataSetTime.timeEntry.time = (int)round(coincTimePS / 1e9);
  //   m_LMFile->write(reinterpret_cast<char*>(&dataSetTime), sizeof(dataSetTime));
  // }
  
  // if(m_oldDetectorPosition != hits.at(0)->getDetectorPosition()) {
  //   m_totalDifferentDetectorPositions++;
  //   std::cout << "new detector number " << m_oldDetectorPosition << " " << hits.at(0)->getDetectorPosition() << std::endl;
  //   ListModeData dataSetPosition;
  //   dataSetPosition.positionEntry.type = 2;
  //   dataSetPosition.positionEntry.oldPosition = m_oldDetectorPosition;
  //   dataSetPosition.positionEntry.newPosition = hits.at(0)->getDetectorPosition();
  //   m_LMFile->write(reinterpret_cast<char*>(&dataSetPosition), sizeof(dataSetPosition));
  //   m_oldDetectorPosition = hits.at(0)->getDetectorPosition();
  // }
  
  //store the normal coincidence
  ListModeData dataSetEvent;
  dataSetEvent.eventEntry.type = 0;
  
  dataSetEvent.eventEntry.ring0 = 0;//hits.at(0)->getRingNumber();
  dataSetEvent.eventEntry.crystal0 = crys0;//hits.at(0)->getCrystalNumber();
  dataSetEvent.eventEntry.layer0 = 0;
  
  dataSetEvent.eventEntry.ring1 = 0;//hits.at(1)->getRingNumber();
  dataSetEvent.eventEntry.crystal1 = crys1;//hits.at(1)->getCrystalNumber();
  dataSetEvent.eventEntry.layer1 = 0;
  
  // const int64_t tDiff = (hits.at(0)->getTime() - hits.at(1)->getTime());
  // if(abs(tDiff) >= 2048)
  //   std::cerr << "CoincidenceToLMTOFExport::processCoincidence; tDiff to large... " << tDiff << std::endl;
  // const uint64_t unTDiff = tDiff + 2048;
  // dataSetEvent.eventEntry.t0MinusT1 = unTDiff & 0xFFF; //to save the signed value in the unsigned field
  
  m_LMFile->write(reinterpret_cast<char*>(&dataSetEvent), sizeof(dataSetEvent));
  if(m_LMFile->fail()) {
    std::cerr << "DiscreteEarth::ToLMFExport; error while writing to file" << std::endl;
    exit(-1);
  }

}

void DiscreteEarth::CloseLMF(){
  
  // close the file
  m_LMFile->close();
  
}
