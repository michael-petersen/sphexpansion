/*
sphmodel.h

functions to read in the spherical model files

MSP 22 Apr 2020 clean version
MSP  3 Jan 2022 make spline flag
MSP  8 Apr 2022 add linear interpolation option


- still needs a linear interpolation fallback

 */
#ifndef SPHMODEL_H
#define SPHMODEL_H

//using namespace std;
using std::cout, std::cerr, std::endl, std::setw, std::vector, std::ifstream, std::ios, std::string, std::ofstream, std::istringstream;

struct SphModel
{
  int NUMR;          // number of entries in the halo table
  vector<double> r;  // the radius array, len NUMR
  vector<double> d;  // the density array, len NUMR
  vector<double> m;  // the enclosed mass array, len NUMR
  vector<double> p;  // the potential array, len NUMR

};


void read_model (string& model_file, SphModel& modeltable) {

  ifstream infile;
  infile.open(model_file);
  if (!infile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
  }

  string line;

  // is numr a dummy variable?
  int numr;
  double rtmp,dtmp,mtmp,ptmp;

  int linenum = 0;

  while (getline(infile, line)) {
    istringstream ss(line);

    if(line.at(0) == '!') {
      // just skip past this
      numr = 0;
    } else {

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
  }

  infile.close();

}


double get_pot (SphModel& sphmodel, double r) {

  // find the closest r value starting from the smallest
  double indx=0;
  while (sphmodel.r[indx] < r) {
    indx++;
  }

  // revert one index: this is the lower bound
  indx--;

  // do interpolation
  double pslope,drt,drp;
  pslope = sphmodel.p[indx+1] - sphmodel.p[indx];
  drt    = sphmodel.r[indx+1] - sphmodel.r[indx];
  drp    = r - sphmodel.r[indx];

  return sphmodel.p[indx] + (pslope/drt)*(drp);

}

double get_density(SphModel& sphmodel, double r) {

  // find the closest r value starting from the smallest
  double indx=0;
  while (sphmodel.r[indx] < r) {
    indx++;
  }

  // revert one index: this is the lower bound
  indx--;

  // do interpolation
  double dslope,drt,drp;
  dslope = sphmodel.d[indx+1] - sphmodel.d[indx];
  drt    = sphmodel.r[indx+1] - sphmodel.r[indx];
  drp    = r - sphmodel.r[indx];

  return sphmodel.d[indx] + (dslope/drt)*(drp);

}

double get_mass (SphModel& sphmodel, double r) {

  // find the closest r value starting from the smallest
  double indx=0;
  while (sphmodel.r[indx] < r) {
    indx++;
  }

  // revert one index: this is the lower bound
  indx--;

  // do interpolation
  double mslope,drt,drp;
  mslope = sphmodel.m[indx+1] - sphmodel.m[indx];
  drt    = sphmodel.r[indx+1] - sphmodel.r[indx];
  drp    = r - sphmodel.r[indx];

  return sphmodel.m[indx] + (mslope/drt)*(drp);

}


#endif
