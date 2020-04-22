/*
sphmodel.h

functions to read in the spherical model files

clean version, MSP 22 April 2020

 */


using namespace std;

#include "spline.h"

struct SphModel
{
  int NUMR;       // number of entries in the halo table
  vector<double> r;  // the radius array, len NUMR
  vector<double> d;  // the density array, len NUMR
  vector<double> m;  // the enclosed mass array, len NUMR
  vector<double> p;  // the potential array, len NUMR

  tk::spline pspline; // spline representation of potential
  tk::spline dspline; // spline representation of density
  tk::spline mspline; // spline representation of mass
};


void read_model (string& model_file, SphModel& modeltable) {

  ifstream infile;
  infile.open(model_file);
  if (!infile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
  }
  
  string line;

  int numr;
  double rtmp,dtmp,mtmp,ptmp;
  
  int linenum = 0;

  while (getline(infile, line)) {
    istringstream ss(line);

    if (linenum==0) {
	ss >> modeltable.NUMR;
	modeltable.r.resize(modeltable.NUMR);
	modeltable.d.resize(modeltable.NUMR);
	modeltable.m.resize(modeltable.NUMR);
	modeltable.p.resize(modeltable.NUMR);
      } else {
        ss >> 
        modeltable.r[linenum-1] >> 
        modeltable.d[linenum-1] >>
        modeltable.m[linenum-1] >>
        modeltable.p[linenum-1];
    }
  linenum ++;

  }

  infile.close();

  // construct spline interpolations
  modeltable.pspline.set_points(modeltable.r,modeltable.p);
  modeltable.dspline.set_points(modeltable.r,modeltable.d);
  modeltable.mspline.set_points(modeltable.r,modeltable.m);

}



// evaluation calls
double get_pot (SphModel& sphmodel, double r) {
  return sphmodel.pspline(r);
}

double get_mass (SphModel& sphmodel, double r) {
  return sphmodel.mspline(r); 
}

double get_density (SphModel& sphmodel, double r) {
  return sphmodel.dspline(r);
}

