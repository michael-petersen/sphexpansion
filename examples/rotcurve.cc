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

// expansion includes
#include "cylexpansion.h"


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

    S->return_forces(coefs,xin, 0., 0.,fx, fy, fz, harmonicflag);

    S->return_density(coefs,xin, 0., 0.,d, harmonicflag);


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

  cout << "Initialising MW disc ... " << endl;
  CylExpansion* MWD = new CylExpansion(cyl_cache_name_mw, cyl_coef_name_mw, cyl_orient_name_mw);

  bool onlymonopole = false;
  MatrixXd mwcoefs,lmccoefs;
  MW->select_coefficient_time(0.0, mwcoefs);
  LMC->select_coefficient_time(0.0, lmccoefs);

  MatrixXd mwdcoscoefs,mwdsincoefs;
  MWD->select_coefficient_time(0.0, mwdcoscoefs, mwdsincoefs);

  double p0,p1,pr,pt,pp;
  MW->determine_fields_at_point_sph(mwcoefs, 0.1,1.570796,0.01,p0,p1,pr,pt,pp);//,0);
  cout << "p0=" << setw(14) << p1 << endl;

  // try returning the whole array
  MatrixXd p1m,prm,ptm,ppm;
  std::tuple<MatrixXd,MatrixXd,MatrixXd,MatrixXd>  X;

  X = MW->determine_weights_at_point_sph(mwcoefs, 0.1,1.570796,0.01);
  p1m = std::get<0>(X);
  prm = std::get<1>(X);
  ptm = std::get<2>(X);
  ppm = std::get<3>(X);
  //cout << "coef0=" << setw(14) << mwcoefs(0,0) << " weight0=" << setw(14) << p1m(0,0) << endl;

  double p1tot=0.0;
  for (int z=0;z<mwcoefs.rows();z++) {
    for (int n=0;n<mwcoefs.cols();n++) {
      p1tot += mwcoefs(z,n)*p1m(z,n);
    }
  }

  cout << "total0=" << setw(14) << p1tot << endl;


  MWD->determine_fields_at_point_cyl(mwdcoscoefs, mwdsincoefs, 0.1,1.570796,0.01,p0,p1,pr,pt,pp,0);//,0);
  cout << "D p0=" << setw(14) << p1 << endl;

  MatrixXd p1mc,prmc,ptmc,ppmc,p1ms,prms,ptms,ppms;
  std::tuple<MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd>  Y;

  Y = MWD->determine_weights_at_point_cyl(mwdcoscoefs, mwdsincoefs, 0.1,1.570796,0.01);
  p1mc = std::get<0>(Y);
  prmc = std::get<1>(Y);
  ptmc = std::get<2>(Y);
  ppmc = std::get<3>(Y);
  p1ms = std::get<4>(Y);
  prms = std::get<5>(Y);
  ptms = std::get<6>(Y);
  ppms = std::get<7>(Y);
  cout << "D coef0=" << setw(14) << mwdcoscoefs(0,0) << " weight0=" << setw(14) << p1mc(0,0) << endl;

  double p1totc=0.0;
  int z = 0;
  //for (int z=0;z<mwcoefs.rows();z++) {
    for (int n=0;n<mwcoefs.cols();n++) {
      p1totc += mwdcoscoefs(z,n)*p1mc(z,n);
    }
  //}

  cout << "D total0=" << setw(14) << p1totc << endl;



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
