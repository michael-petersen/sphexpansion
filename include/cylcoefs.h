/*
cylecoefs.h

functions to handle the preparations for the cylinder coefficient files

MSP 5 May 2020 clean version

now takes raw outcoef files -- need helper functions to check for data corruption from restarts, etc.

todo:
  - convert to eigen arrays

 */
#ifndef CYLCOEFS_H
#define CYLCOEFS_H


using namespace std;

// create 2- and 3-d array types
typedef boost::multi_array<double, 3> array_type3;
typedef boost::multi_array<double, 2> array_type2;

struct CylCoefs
{
  int MMAX;            // the number of azimuthal harmonics
  int NORDER;            // the number of radial terms
  int NUMT;            // the number of timesteps
  
  vector<double> t;    // the time, len NUMT  

  array_type3 coscoefs;   // the cosine coefficient table, sized NUMT,MMAX+1,NORDER
  array_type3 sincoefs;   // the   sine coefficient table, sized NUMT,MMAX+1,NORDER

};


void select_coefficient_time(double desired_time,
			     CylCoefs coeftable,
			     array_type2& coscoefs_at_time,
			     array_type2& sincoefs_at_time) {
  /*
    linear interpolation to get the coefficient matrix at a specific time

   time units must be virial time units to match the input coefficient table
   */

  // coeftable.t is assumed to be evenly spaced
  double dt = coeftable.t[1] - coeftable.t[0];

  int indx = (int)( (desired_time-coeftable.t[0])/dt);

  // guard against wanton extrapolation: should this stop the model?
  if (indx<0) cerr << "select_coefficient_time: time prior to simulation start selected. setting to earliest step." << endl;
  if (indx>coeftable.NUMT-2) cerr << "select_coefficient_time: time after to simulation end selected. setting to latest step." << endl;

  if (indx<0) indx = 0;
  if (indx>coeftable.NUMT-2) indx = coeftable.NUMT - 2;

  double x1 = (coeftable.t[indx+1] - desired_time)/dt;
  double x2 = (desired_time - coeftable.t[indx])/dt;

  coscoefs_at_time.resize(boost::extents[coeftable.MMAX+1][coeftable.NORDER]);
  sincoefs_at_time.resize(boost::extents[coeftable.MMAX+1][coeftable.NORDER]);

  for (int m=0; m<=coeftable.MMAX; m++){
    for (int n=0; n<coeftable.NORDER; n++) {
      coscoefs_at_time[m][n] = (x1*coeftable.coscoefs[indx][m][n] + x2*coeftable.coscoefs[indx+1][m][n]);

      if (m) sincoefs_at_time[m][n] = (x1*coeftable.sincoefs[indx][m][n] + x2*coeftable.sincoefs[indx+1][m][n]);
    }
  }
  
}


void read_coef_file (string& coef_file, CylCoefs& coeftable) {

  ifstream in(coef_file.c_str());

  // read a template version first, then reread with NUMT specified

  double tnow;
  // first thing in is NUMT,LMAX,NORDER
  in.read((char *)&tnow, sizeof(double));
  in.read((char *)&coeftable.MMAX, sizeof(int));
  in.read((char *)&coeftable.NORDER, sizeof(int));

  int end,tmp;
  in.seekg (0, ios::end);
  end = in.tellg();

  // compute the number of full timesteps: extra 8s are for tnow,MMAX,NORDER ahead of every coefficient set
  coeftable.NUMT = end/(((coeftable.MMAX+coeftable.MMAX+1)*coeftable.NORDER)*sizeof(double)  + 8 + 8);

  //cout << setw(14) << (coeftable.MMAX*coeftable.NORDER) << setw(14) << sizeof(double) << setw(14) << end << endl;
  //cout << coeftable.NUMT << endl;

  cout << "cylcoefs.read_coef_file: reading NUMT, LMAX, NORDER from file . . . " << endl;
  cout << setw(18) << coeftable.NUMT << setw(18) << coeftable.MMAX <<
    setw(18) << coeftable.NORDER << endl;

  // reset to the beginning
  in.seekg (0, ios::beg);


  // resize the coefs array appropriately
  coeftable.coscoefs.resize(boost::extents[coeftable.NUMT][coeftable.MMAX+1][coeftable.NORDER]);
  coeftable.sincoefs.resize(boost::extents[coeftable.NUMT][coeftable.MMAX+1][coeftable.NORDER]);
  coeftable.t.resize(coeftable.NUMT);
  
  // now cycle through each time
  for (int tt=0;tt<coeftable.NUMT;tt++) {
  
    in.read((char *)&coeftable.t[tt], sizeof(double));
    in.read((char *)&tmp, sizeof(int));
    in.read((char *)&tmp, sizeof(int));

    // debug outcoef file
    //cout << "tnow=" << coeftable.t[tt] << " norder=" << tmp << endl;

    for (int m=0; m<=coeftable.MMAX; m++) {
      
      for (int ir=0; ir<coeftable.NORDER; ir++) {
        in.read((char *)&coeftable.coscoefs[tt][m][ir], sizeof(double));
      }

      if (m) {
	for (int ir=0; ir<coeftable.NORDER; ir++) {
          in.read((char *)&coeftable.coscoefs[tt][m][ir], sizeof(double));
        }
      }
      
    } // MMAX loop
  } // NUMT loop



  cout << "success!!" << endl;

}



#endif
