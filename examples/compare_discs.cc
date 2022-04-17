/*
example code:
compare two different disc cache files

compile string:
clang++ --std=c++17 -lyaml-cpp -I/opt/local/include  -I../include/ fullintegrate.cc -o obj/fullintegrate

MSP 17 Apr 2022 initial commit

*/

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

// the cylindrical expansion headers
#include "cylexpansion.h"



int main () {

  // MW with high-resolution basis
  cout << "Initialising MW disc ... " << endl;
  CylExpansion* MWD = new CylExpansion(cyl_cache_name_mw, cyl_coef_name_mw, cyl_orient_name_mw);

  // get coefficients...shared between both cache files
  MatrixXd mwdcoscoefs,mwdsincoefs;
  MWD->select_coefficient_time(0.0, mwdcoscoefs, mwdsincoefs);

  // MW with high-resolution basis
  cout << "Initialising HI-RES MW disc ... " << endl;
  string cyl_cache_name_mw_hires = "../data/disc.cache.run9mlde";
  CylExpansion* MWDH = new CylExpansion(cyl_cache_name_mw_hires, cyl_coef_name_mw, cyl_orient_name_mw);




  double r2tmp = 0.004;
  double phitmp = 0.0;
  double zvir = -0.001;

  CylForce ftable1,ftable2;
  MWD->get_table_forces(r2tmp,zvir,ftable1);
  MWDH->get_table_forces(r2tmp,zvir,ftable2);

  cout << "Force table check:" << setw(14) << ftable1.rforceC(0,0) << setw(14) << ftable2.rforceC(0,0) << std::endl;

  double tpotl0,tpotl,fr,fp,fztmp;
  MWD->determine_fields_at_point_cyl(mwdcoscoefs,mwdsincoefs,
				     r2tmp,phitmp,zvir,
				     tpotl0,tpotl,
				     fr,fp,fztmp,false,false,false);

  cout << "Forces (lowres):" << setw(14) << fr << setw(14) << fp << setw(14) << fztmp << endl;

  MWDH->determine_fields_at_point_cyl(mwdcoscoefs,mwdsincoefs,
             r2tmp,phitmp,zvir,
             tpotl0,tpotl,
             fr,fp,fztmp,false,false,false);

  cout << "Forces ( hires):" << setw(14) << fr << setw(14) << fp << setw(14) << fztmp << endl;


}
