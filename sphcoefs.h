/*
sphcoefs.h

functions to handle the preparations for the coefficient files

clean version, MSP 22 April 2020

 */

using namespace std;

// create 2- and 3-d array types
typedef boost::multi_array<double, 3> array_type3;
typedef boost::multi_array<double, 2> array_type2;


struct SphCoefs
{
  int LMAX;            // the number of azimuthal harmonics
  int NMAX;            // the number of radial terms
  int NUMT;            // the number of timesteps
  
  vector<double> t;    // the time, len NUMT  

  array_type3 coefs;   // the coefficient table, sized NUMT,(LMAX+1)*(LMAX+1),NMAX

};
 


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

  cout << "success!!" << endl;

}
