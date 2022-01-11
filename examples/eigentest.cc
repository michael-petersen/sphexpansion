/*
tests for the initial MW using eigen

compile string:
clang++ --std=c++17 -I/opt/local/include -L/opt/local/lib -I../include/ -I/opt/local/include/eigen3 eigentest.cc -o obj/eigentest

MSP 22 Dec 2021 first version

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdio.h>

// override eigen flag for this test
#undef HAVEEIGEN
#define HAVEEIGEN 1

// eigen includes
#include <Eigen/Dense>

// expansion includes
#include "expansion.h"



void make_rotation_curve(SphExpansion* S,
			 array_type2 coefs,
			 double xmin,
			 double xmax,
			 int nsamples,
			 string outfile,
			 bool monopole=false)
{
  /*
  xmin and xmax come in PHYSICAL units
  and the rotation curve is returned in physical units
   */

  if (monopole) cout << "Only using monopole for rotation curve." << endl;

  ofstream mwrotation;
  mwrotation.open(outfile);

  mwrotation << "# radius [kpc]; vcirc [km/s] ; f_x [km/s/s] ; f_y [km/s/s] ; f_z [km/s/s];" << endl;

  double dx = (xmax-xmin)/nsamples;
  double xin;
  double fx,fy,fz,d,physdens;

  // demo to make a rotation curve
  for (int xx=0; xx<nsamples; xx++) {

    // the location in inertial space of the points to check (yin=zin=0, just checking x-axis right now)
    xin = xx*dx + xmin; // in kpc

    S->return_forces(
		  coefs,
		  xin, 0., 0.,
		     fx, fy, fz, monopole);

    S->return_density(
		  coefs,
		  xin, 0., 0.,
		     d, monopole);


    virial_to_physical_density(d, physdens);

    mwrotation << setw(14) << xin << setw(14) << sqrt(xin*-fx) <<
                  setw(14) << fx << setw(14) << fy << setw(14) << fz <<
                  setw(14) << d << setw(14) << physdens << endl;

  }

  mwrotation.close();

}



int main () {

  // MW
  cout << "Initialising MW ... " << endl;
  string sph_cache_name_mw  =     "../data/run068s22h/SLGridSph.cache.mw.run068s22h";
  string model_file_mw      =     "../data/run068s22h/SLGridSph.mw.s22h";
  string coef_file_mw       =     "../data/run068s22h/simpleoutcoef.nofac.mw.run068s22h";
  string orient_file_mw     =     "../data/run068s22h/mw.orient.run068s22h.smth";

  SphExpansion* MW;
  MW = new SphExpansion(sph_cache_name_mw, model_file_mw, coef_file_mw, orient_file_mw);


  bool onlymonopole = false;
  array_type2 mwcoefs;
  MW->select_coefficient_time(0.0, mwcoefs);

  int numW=100;

  Eigen::MatrixXd self_grav_coefs;
  //MW->get_selfgravity_coefficients_eigen(self_grav_coefs);

  Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(numW, numW);

  ret.resize(40,40);

	Eigen::MatrixXd factrl;
  factorial_eigen(6, factrl);


  string rotationfile="tests/MWrotation.txt";

  make_rotation_curve(MW,
		      mwcoefs,
		      0.1,
		      120.,
		      50,
		      rotationfile, true);


}
