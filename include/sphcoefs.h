/*
sphcoefs.h

functions to handle the preparations for the coefficient files

MSP 22 Apr 2020 clean version
MSP 23 Apr 2020 add coefficient interpolation
MSP 15 Sep 2020 improve debug, fix pre-simulation interpolation
MSP 14 May 2021 added self-gravity normalisation coefficients
MSP 24 Dec 2021 change some debug outputs and spline calls to preprocessor flags
MSP  9 Apr 2022 converted to Eigen


 */
#ifndef SPHCOEFS_H
#define SPHCOEFS_H

#include <Eigen/StdVector>
#include <Eigen/Dense>


// standard namespace includes
using std::cout, std::cerr, std::endl, std::setw, std::vector, std::ifstream, std::ios, std::string, std::ofstream, std::istringstream;

// Eigen MatrixXd
using Eigen::MatrixXd;


struct SphCoefs
{
  int LMAX;            // the number of azimuthal harmonics
  int NMAX;            // the number of radial terms
  int NUMT;            // the number of timesteps

  vector<double> t;    // the time, len NUMT

  std::vector<MatrixXd> coefs;   // the coefficient table, sized NUMT,(LMAX+1)*(LMAX+1),NMAX

};


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

  // resize the vector for time
  coeftable.coefs.resize(coeftable.NUMT);
  coeftable.t.resize(coeftable.NUMT);

  // now cycle through each time
  for (int tt=0;tt<coeftable.NUMT;tt++) {

    // resize the matrix
    coeftable.coefs[tt].resize(numl,coeftable.NMAX);

    in.read((char *)&coeftable.t[tt], sizeof(double));

    for (int l=0; l<numl; l++) {
      for (int ir=0; ir<coeftable.NMAX; ir++) {
        in.read((char *)&coeftable.coefs[tt](l,ir), sizeof(double));
      }
    }
  }


  cout << "success!!" << "\n";

}



#endif
