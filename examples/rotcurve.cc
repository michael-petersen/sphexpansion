/*
tests for the initial MW and LMC
- generate a rotation curve
- integrate a circular orbit

compile string:
clang++ --std=c++17 -I/opt/local/include -L/opt/local/lib -I../include/ rotcurve.cc -o obj/rotcurve

MSP 24 Apr 2020 first version
MSP 29 Sep 2020 use as a test for density
MSP 28 Sep 2021 fix for modern compatibility
MSP 11 Jan 2022 convert to c++17 standard
MSP 16 Apr 2022 point to stable example files

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

// expansion includes
#include "sphexpansion.h"



void make_rotation_curve(SphExpansion* S,
			 MatrixXd coefs,
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
	int harmonicflag = 2047; // default; up to l<=10
  if (monopole) {
		std::cout << "Only using monopole for rotation curve." << endl;
		harmonicflag = 0;
	}

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
		     fx, fy, fz, false,false,false,1000, harmonicflag);

    S->return_density(
		  coefs,
		  xin, 0., 0.,
		     d, false,false,false,1000, harmonicflag);


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
  SphExpansion* MW = new SphExpansion(sph_cache_name_mw, model_file_mw, coef_file_mw, orient_file_mw);

  // LMC
  cout << "Initialising LMC ... " << endl;
  SphExpansion* LMC = new SphExpansion(sph_cache_name_lmc, model_file_lmc, coef_file_lmc, orient_file_lmc);


  bool onlymonopole = false;
  MatrixXd mwcoefs,lmccoefs;
  MW->select_coefficient_time(0.0, mwcoefs);
  LMC->select_coefficient_time(0.0, lmccoefs);

  vector<double> velcentre(3),centre(3);
  return_vel_centre(0.0, LMC->orient, velcentre);
  cout << "C.O.M. velocity vx=" << velcentre[0] << " vy=" <<
    velcentre[1] << " vz=" << velcentre[2] << endl;

  return_centre(0.0, LMC->orient, centre);
  cout << "C.O.M. x=" << centre[0] << " y=" <<
    centre[1] << " z=" << centre[2] << endl;

  return_vel_centre(-1.0, LMC->orient, velcentre);
  cout << "C.O.M. velocity vx=" << velcentre[0] << " vy=" <<
    velcentre[1] << " vz=" << velcentre[2] << endl;

  return_centre(-1.0, LMC->orient, centre);
  cout << "C.O.M. x=" << centre[0] << " y=" <<
    centre[1] << " z=" << centre[2] << endl;

  return_vel_centre(-2.0, LMC->orient, velcentre);
  cout << "C.O.M. velocity vx=" << velcentre[0] << " vy=" <<
    velcentre[1] << " vz=" << velcentre[2] << endl;

  return_centre(-2.0, LMC->orient, centre);
  cout << "C.O.M. x=" << centre[0] << " y=" <<
    centre[1] << " z=" << centre[2] << endl;



  string rotationfile="tests/MWrotation.txt";

  make_rotation_curve(MW,
		      mwcoefs,
		      0.1,
		      120.,
		      50,
		      rotationfile, true);

  rotationfile="tests/LMCrotation.txt";

  make_rotation_curve(LMC,
		      lmccoefs,
		      0.1,
		      120.,
		      50,
		      rotationfile, true);

}
