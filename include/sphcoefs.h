/*
sphcoefs.h

functions to handle the preparations for the coefficient files

MSP 22 Apr 2020 clean version
MSP 23 Apr 2020 add coefficient interpolation 
MSP 15 Sep 2020 improve debug, fix pre-simulation interpolation

notes
-the answer is yes, spline is expensive. moving the call outside helps; a more simple interpolation may help even further.

 */

#define DEEPDEBUGCOEFS 0

using namespace std;

// set up for spline, specify if needed
#include "spline.h"
bool splinecoefs = false;

bool debugcoefs = false;


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

void select_coefficient_time(double desired_time, SphCoefs coeftable,
			     array_type2& coefs_at_time, int order=-1);


void select_coefficient_time(double desired_time, SphCoefs coeftable,
			     array_type2& coefs_at_time, int order) {
  /*
    linear interpolation to get the coefficient matrix at a specific time

   time units must be virial time units to match the input coefficient
   table

   if order<0, will select all orders. if order>0, will select only
   specified order.
   */

  int numl;

  numl = (coeftable.LMAX+1)*(coeftable.LMAX+1);

  coefs_at_time.resize(boost::extents[numl][coeftable.NMAX]);  
  
  // coeftable.t is assumed to be evenly spaced
  double dt = coeftable.t[1] - coeftable.t[0];

  int indx = (int)( (desired_time-coeftable.t[0])/dt);

  // guard against wanton extrapolation: should this stop the model?
  //if (indx<0) cerr << "select_coefficient_time: time prior to simulation start selected. setting to earliest step." << endl;

  // guard against going past the end of the simulation
  if (indx>coeftable.NUMT-2) cerr << "select_coefficient_time: time after to simulation end selected. setting to latest step." << endl;

#if DEEPDEBUGCOEFS
  cout << indx << endl;
#endif
  
  if (indx<0) {


    for (int l=0; l<numl; l++){
      for (int n=0; n<coeftable.NMAX; n++) {
        coefs_at_time[l][n] = coeftable.coefs[0][l][n];
      }
    }
    
  } else {

    // case where the simulations are not before the beginning of the simulation
    // interpolate from the two closest times

    double x1 = (coeftable.t[indx+1] - desired_time)/dt;
    double x2 = (desired_time - coeftable.t[indx])/dt;

#if DEEPDEBUGCOEFS
    cout << "dt=" << setw(16) << dt << "  t[indx+1]="  << setw(16) << coeftable.t[indx+1]
	                            << "  t[indx]="  << setw(16) << coeftable.t[indx]
                                    << "  desiredT="  << setw(16) << desired_time << endl;
    cout << "x1/x2=" << setw(16) << x1 << setw(14) << x2 << endl;
#endif

    for (int l=0; l<numl; l++){
      for (int n=0; n<coeftable.NMAX; n++) {
        coefs_at_time[l][n] = (x1*coeftable.coefs[indx][l][n] + x2*coeftable.coefs[indx+1][l][n]);
      }
    }
  }

  // go through and zero out non-selected orders

  if (order>0) {
    for (int l=0; l<numl; l++){
    
      if (l == order) continue;
      
      for (int n=0; n<coeftable.NMAX; n++) {
        coefs_at_time[l][n] = 0.0;
      }
      
    }
  }
  
}




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

  cout << "sphcoefs.read_coef_file: reading NUMT, LMAX, NMAX from file . . . " << endl;
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

  if (splinecoefs) {
  cout << "setting up coefficient interpolation . . . ";
  make_coef_splines(coeftable);
  }

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

  cout << "sphcoefs.read_coef_file: reading coefficients from file . . . ";

  if (debugcoefs) {
  cout << endl << "sphcoefs.read_coef_file: reading NUMT, LMAX, NMAX from file . . . " << endl;
  cout << setw(18) << coeftable.NUMT << setw(18) << coeftable.LMAX <<
    setw(18) << coeftable.NMAX << endl;
  }

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

  if (splinecoefs) {
  cout << "setting up coefficient interpolation . . . ";
  make_coef_splines(coeftable);
  }

  cout << "success!!" << endl;

}



