/*
sphcoefs.h

functions to handle the preparations for the coefficient files

MSP 22 Apr 2020 clean version
MSP 23 Apr 2020 add coefficient interpolation
MSP 15 Sep 2020 improve debug, fix pre-simulation interpolation
MSP 14 May 2021 added self-gravity normalisation coefficients
MSP 24 Dec 2021 change some debug outputs and spline calls to preprocessor flags


notes
-the answer is yes, spline is expensive. moving the call outside helps; a more simple interpolation may help even further.

 */
#ifndef SPHCOEFS_H
#define SPHCOEFS_H


//#include "basis.h"

#if HAVEEIGEN
#include <Eigen/Dense>
#endif


//using namespace std;
using std::cout, std::cerr, std::endl, std::setw, std::vector, std::ifstream, std::ios, std::string, std::ofstream, std::istringstream;


// create 2- and 3-d array types
typedef boost::multi_array<double, 3> array_type3;
typedef boost::multi_array<double, 2> array_type2;

#if SPLINECOEFS
// create a spline array for the coefficients
// set up for spline, specify if needed
#include "spline.h"
typedef boost::multi_array<tk::spline, 2> spline_array;
#endif

struct SphCoefs
{
  int LMAX;            // the number of azimuthal harmonics
  int NMAX;            // the number of radial terms
  int NUMT;            // the number of timesteps

  vector<double> t;    // the time, len NUMT

  array_type3 coefs;   // the coefficient table, sized NUMT,(LMAX+1)*(LMAX+1),NMAX

#if SPLINECOEFS
  spline_array coefsplines; // spline version of the coefficients, sized (LMAX+1)*(LMAX+1),NMAX
#endif

};


#if SPLINECOEFS
void spline_coefficient_time(double desired_time, SphCoefs coeftable, array_type2& coefs_at_time) {
  /*
    interpolate to get the coefficient matrix at a specific time
   */

  int numl = (coeftable.LMAX+1)*(coeftable.LMAX+1);

  coefs_at_time.resize(boost::extents[numl][coeftable.NMAX]);

  for (int l=0; l<numl; l++){

    for (int n=0; n<coeftable.NMAX; n++) {

      coefs_at_time[l][n] = coeftable.coefsplines[l][n](desired_time);

    }

  }

}




void make_coef_splines( SphCoefs& coeftable) {
  /*
  step through the coefficients and make splines.
   */

  // initialise an empty array to hold the coefficients per order
  vector<double> tmparray(coeftable.NUMT);

  int numl = (coeftable.LMAX+1)*(coeftable.LMAX+1);

  coeftable.coefsplines.resize(boost::extents[numl][coeftable.NMAX]);

  for (int l=0; l<numl; l++){

    for (int n=0; n<coeftable.NMAX; n++) {

      // recast to a vector: this may be a spot to speed up computation if needed?
      for (int t=0; t<coeftable.NUMT; t++) tmparray[t] = coeftable.coefs[t][l][n];

      // construct the splines
      coeftable.coefsplines[l][n].set_points(coeftable.t,tmparray);

    }

  }

}
#endif


void read_coef_file_raw(string& coef_file, SphCoefs& coeftable) {

  /*
    read in the self-describing coefficient file

  */
  ifstream in(coef_file.c_str());

  double tnow,scale;

  char buf[64];
  in.read((char *)&buf, 64*sizeof(char));
  in.read((char *)&tnow, sizeof(double));
  in.read((char *)&scale, sizeof(double));
  in.read((char *)&coeftable.NMAX, sizeof(int));
  in.read((char *)&coeftable.LMAX, sizeof(int));

  double tmp;
  int end,tmpint;
  in.seekg (0, ios::end);
  end = in.tellg();

  int numl = (coeftable.LMAX+1) * (coeftable.LMAX+1);

  // compute the number of full timesteps: extra 8s are for tnow,MMAX,NORDER ahead of every coefficient set
  coeftable.NUMT = end/((numl*coeftable.NMAX)*sizeof(double)  + 8 + 16 + 64);

  cout << "sphcoefs::read_coef_file_raw: reading NUMT, LMAX, NMAX from file . . . " << endl;
  cout << setw(18) << coeftable.NUMT << setw(18) << coeftable.LMAX <<
    setw(18) << coeftable.NMAX << endl;

  // resize the coefs array appropriately
  coeftable.coefs.resize(boost::extents[coeftable.NUMT][numl][coeftable.NMAX]);
  coeftable.t.resize(coeftable.NUMT);

  // now cycle through each time
  for (int tt=0;tt<coeftable.NUMT;tt++) {

    in.read((char *)&buf, 64*sizeof(char));
    in.read((char *)&coeftable.t[tt], sizeof(double));
    in.read((char *)&tmp, sizeof(double));
    in.read((char *)&tmpint, sizeof(int));
    in.read((char *)&tmpint, sizeof(int));

    for (int l=0; l<numl; l++) {
      for (int ir=0; ir<coeftable.NMAX; ir++) {
        in.read((char *)&coeftable.coefs[tt][l][ir], sizeof(double));
      }
    }
  }

#if SPLINECOEFS
  cout << "sphcoefs::read_coef_file_raw: setting up coefficient interpolation . . . ";
  make_coef_splines(coeftable);
#endif

  cout << "success!!" << endl;

}



void read_coef_file (string& coef_file, SphCoefs& coeftable) {

  /*
    read in the self-describing coefficient file

  */
  ifstream in(coef_file.c_str());
  if (!in) {
        cout << "sphcoefs::read_coef_file: Unable to open file!\n";
        exit(1);
  }

  // first thing in is NUMT,LMAX,NMAX
  in.read((char *)&coeftable.NUMT, sizeof(int));
  in.read((char *)&coeftable.LMAX, sizeof(int));
  in.read((char *)&coeftable.NMAX, sizeof(int));

  std::cout << "sphcoefs::read_coef_file: reading coefficients from file . . . ";

#if DEBUGCOEFS
  std::cout << "\n" << "sphcoefs::read_coef_file: reading NUMT, LMAX, NMAX from file . . . " << "\n";
  std::cout << setw(18) << coeftable.NUMT << setw(18) << coeftable.LMAX << setw(18) << coeftable.NMAX << "\n";
#endif

  // resize the coefs array appropriately
  int numl = (coeftable.LMAX+1) * (coeftable.LMAX+1);
  coeftable.coefs.resize(boost::extents[coeftable.NUMT][numl][coeftable.NMAX]);
  coeftable.t.resize(coeftable.NUMT);

  // now cycle through each time
  for (int tt=0;tt<coeftable.NUMT;tt++) {

    in.read((char *)&coeftable.t[tt], sizeof(double));

    for (int l=0; l<numl; l++) {
      for (int ir=0; ir<coeftable.NMAX; ir++) {
        in.read((char *)&coeftable.coefs[tt][l][ir], sizeof(double));
      }
    }
  }

#if SPLINECOEFS
  cout << "sphcoefs::read_coef_file: setting up coefficient interpolation . . . ";
  make_coef_splines(coeftable);
#endif

  cout << "success!!" << "\n";

}



void get_selfgravity_coefficients(SphCoefs coeftable,
			          array_type3& self_grav_coefs, int lorder=-1, int norder=-1, bool monopolenorm=false)
{
  // function to produce coefficients that have been normalised such that the self-gravity can be compared
  // this means including (1/(4*pi)) * (2*l+1) * ((l-m)!/(l+m)!)
  // this will set up another copy of the coefficients, so be warned: this could be large

  // todo:
  // add options for single lorder,morder return
  // add monopolenorm option

  // the full version is in expansion.h, where the tables are exposed.

  int numl = coeftable.LMAX;
  int numn = coeftable.NMAX;
  int l,loffset,moffset,m,n,t;

  self_grav_coefs.resize(boost::extents[coeftable.NUMT][(coeftable.LMAX+1)*(coeftable.LMAX+1)][coeftable.NMAX]);

  double fac1,fac2;

  fac1 = 0.25/M_PI;



#if HAVEEIGEN
  Eigen::MatrixXd factrl;
  factorial_eigen(numl, factrl);
#else
  array_type2 factrl;
  factorial(numl, factrl);
#endif

  for (t=0;t<coeftable.NUMT;t++) {

    for (n=0;n<numn;n++) self_grav_coefs[t][0][n] = fac1*coeftable.coefs[t][0][n];

    for (l=1, loffset=1; l<=numl; loffset+=(2*l+1), l++) {
      for (m=0, moffset=0; m<=l; m++) {
        fac1 = (2.0*l+1.0)/(4.0*M_PI);
        if (m==0) {
          for (n=0;n<numn;n++) self_grav_coefs[t][loffset+moffset][n] = fac1*coeftable.coefs[t][loffset+moffset][n];
        } else {
#if HAVEEIGEN
	  fac2 = 2.0 * fac1 * factrl(l,m);
#else
    fac2 = 2.0 * fac1 * factrl[l][m];
#endif
          for (n=0;n<numn;n++) self_grav_coefs[t][loffset+moffset][n] = fac2*coeftable.coefs[t][loffset+moffset][n];
        }
      }
    }
  }

}

#endif
