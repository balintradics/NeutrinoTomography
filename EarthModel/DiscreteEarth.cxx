#include "math.h"
#include "iostream"
#include "fstream"

#include "DiscreteEarth.h"

using namespace std;

DiscreteEarth::DiscreteEarth(float dcell){

  m_DCell = dcell;// [km]
  if (m_DCell < 50.0){
    cout << "Cell size is too small, to save memory resizing to 50 km" << endl;
    m_DCell = 50.0;
  }
  
  // How many cells we need to fill the full Earth volume?
  // We calculate it dynamically...
  m_NCells = 0;
  m_EarthCells = NULL;
  for(float ix = -R_E; ix <= R_E; ix = ix + m_DCell){
    for(float iy = -R_E; iy <= R_E; iy = iy + m_DCell){
      for(float iz = -R_E; iz <= R_E; iz = iz + m_DCell){
	// skip cells that would be outside of the Earth volume
	if(sqrt(ix*ix + iy*iy + iz*iz) > R_E){
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
  float rho_ret = 2.6;
  float r = 0;
  float R = 0;
  for(float ix = -R_E; ix <= R_E; ix = ix + m_DCell){
    for(float iy = -R_E; iy <= R_E; iy = iy + m_DCell){
      for(float iz = -R_E; iz <= R_E; iz = iz + m_DCell){
	r = sqrt(ix*ix+iy*iy+iz*iz);
	if(r > R_E){
	}else{
	  Cell_t cell;
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
	  
	  cell.rho = rho_ret;
	  cell.x = ix;
	  cell.y = iy;
	  cell.z = iz;
	  
	  m_EarthCells[Celli] = cell;
	  Celli++;
	}
      }
    }
  }

  cout << "Allocated " << m_NCells << ", " << Celli << " cells" << endl;
  
}

DiscreteEarth::~DiscreteEarth(){

  delete[] m_EarthCells; 

}

void DiscreteEarth::SaveDensityRToCSV(const char * ofilename){
    
  ofstream outfile;
  outfile.open(ofilename);
  float r = 0;
  for(unsigned int i = 0; i < m_NCells; i++){
    // only print along x = 0, y = 0
    if( fabs(m_EarthCells[i].x) <= m_DCell && fabs(m_EarthCells[i].y) <= m_DCell && m_EarthCells[i].z >= 0){
      r = sqrt(m_EarthCells[i].x*m_EarthCells[i].x+m_EarthCells[i].y*m_EarthCells[i].y+m_EarthCells[i].z*m_EarthCells[i].z);
      outfile  << m_EarthCells[i].z << "\t" << m_EarthCells[i].rho << endl;
    }
  }
  outfile.close();

}

void DiscreteEarth::SaveEarthToCSV(const char * ofilename){
  cout << "Saving : " << m_NCells << " cells values to " << ofilename << endl;    
  ofstream outfile;
  outfile.open(ofilename);
  for(unsigned int i = 0; i < m_NCells; i++){
    outfile  << m_EarthCells[i].x << ", " << m_EarthCells[i].y << ", " << m_EarthCells[i].z << ", " << m_EarthCells[i].rho << endl;
  }
  outfile.close();

}
// latitude, longitude, flux
void DiscreteEarth::SaveFluxToCSV(const char * ofilename){
}
