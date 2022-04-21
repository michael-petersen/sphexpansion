/*
header file for fullintegrate.cc capabilities

MSP 21 Apr 2022 cleaned version v0.2.3

*/
#ifndef FULLINTEGRATE_H
#define FULLINTEGRATE_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdio.h>

// eigen includes
#include <Eigen/Dense>
using Eigen::MatrixXd;

// set model parameters
#include "modelfiles.h"

// the spherical expansion headers
#include "sphexpansion.h"

// the cylindrical expansion headers
#include "cylexpansion.h"


class MWLMC
{
private:

  void initialise();

public:

  // the expansion classes: these might be better private
  SphExpansion* MW;
  SphExpansion* LMC;
  CylExpansion* MWD;

  // the constructor (no type, no arguments, all defined in modelfiles.h)
  MWLMC();

  // return all fields for the MW halo (in the frame of the MW halo)
  void mwhalo_fields(MatrixXd mwcoefs,
                     double t, double x, double y, double z,
                     double& fx, double& fy, double& fz, double& pot, double& dens,
                     int mwhharmonicflag=127, bool verbose=false);

  // return all fields for the LMC halo
  void lmc_fields(MatrixXd lmccoefs,
                  double t, double x, double y, double z,
                  double& fx, double& fy, double& fz, double& pot, double& dens,
                  int lmcarmonicflag=127, bool verbose=false);

  // return all fields for the MW disc
  // NOTE: density does not work here. not enabled yet for cylindrical expansions.
  //       leave as an inspirational placeholder
  void mwd_fields(MatrixXd mwdcoscoefs, MatrixXd mwdsincoefs,
                  double t, double x, double y, double z,
                  double& fx, double& fy, double& fz, double& pot, double& dens,
                  int mwdharmonicflag=127,
                  bool verbose=false);

  // return total forces
  // @IMPROVE: write all potential as well
  void all_forces(MatrixXd mwcoefs, MatrixXd lmccoefs, MatrixXd mwdcoscoefs, MatrixXd mwdsincoefs,
                  double t, double x, double y, double z,
                  double& fx, double& fy, double& fz,
                  int mwhharmonicflag=127, int mwdharmonicflag=127, int lmcharmonicflag=127,
                  bool verbose=false);

  // compute an orbit integration in all three components
  void orbit(vector<double> xinit,
             vector<double> vinit,
             int nint,
             double dt,
             MatrixXd& orbit,
             int mwhharmonicflag=127, int mwdharmonicflag=127, int lmcharmonicflag=127,
             bool fixedtime=false);

  // print an orbit array
  void print_orbit(MatrixXd orbit, string orbitfile);
};

#endif
