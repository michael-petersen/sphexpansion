/*
sphcoefs.h

functions to handle the preparations for the coefficient files

clean version, MSP 22 April 2020

add coefficient interpolation MSP 23 April 2020

 */

using namespace std;

#include "spline.h"

// create 2- and 3-d array types
typedef boost::multi_array<double, 3> array_type3;
typedef boost::multi_array<double, 2> array_type2;

// create a spline array for the coefficients
typedef boost::multi_array<tk::spline, 2> spline_array;


struct SphCoefs
{
  int LMAX;            // the number of azimuthal harmonics
  int NMAX;            // the number of radial terms
  int NUMT;            // the number of timesteps
  
  vector<double> t;    // the time, len NUMT  

  array_type3 coefs;   // the coefficient table, sized NUMT,(LMAX+1)*(LMAX+1),NMAX

  spline_array coefsplines; // spline version of the coefficients, sized (LMAX+1)*(LMAX+1),NMAX

};


void select_coefficient_time(double desired_time, SphCoefs coeftable, array_type2& coefs_at_time) {
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


void read_coef_file (string& coef_file, SphCoefs& coeftable) {

  /*
    read in the self-describing coefficient file

  */
  ifstream in(coef_file.c_str());

  // first thing in is NUMT,LMAX,NMAX
  in.read((char *)&coeftable.NUMT, sizeof(int));
  in.read((char *)&coeftable.LMAX, sizeof(int));
  in.read((char *)&coeftable.NMAX, sizeof(int));

  cout << "sphcoefs.read_coef_file: reading NUMT, LMAX, NMAX from file . . . " << endl;
  cout << setw(18) << coeftable.NUMT << setw(18) << coeftable.LMAX <<
    setw(18) << coeftable.NMAX << endl;

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

  cout << "setting up coefficient interpolation . . . ";
  make_coef_splines(coeftable);

  cout << "success!!" << endl;

}



