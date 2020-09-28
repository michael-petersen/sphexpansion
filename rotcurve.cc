/*
tests for the initial MW and LMC
- generate a rotation curve
- integrate a circular orbit

compile string: 
clang++ -I/opt/local/include -L/opt/local/lib -Iinclude/ rotcurve.cc -o rotcurve

MSP 24 Apr 2020 first version

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdio.h>

// boost includes
#include "boost/multi_array.hpp"

// expansion includes
#include "expansion.h"

// integration include
#include "leapfrog.h"

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
  double fx,fy,fz,d;

  // demo to make a rotation curve
  for (int xx=0; xx<nsamples; xx++) {

    // the location in inertial space of the points to check (yin=zin=0, just checking x-axis right now)
    xin = xx*dx + xmin; // in kpc

    S->return_forces(S,
		  coefs,
		  xin, 0., 0.,
		     fx, fy, fz, true);

    S->return_density(S,
		  coefs,
		  xin, 0., 0.,
		     d, true);
    
    mwrotation << setw(14) << xin << setw(14) << sqrt(xin*-fx) <<
                  setw(14) << fx << setw(14) << fy << setw(14) << fz <<
                  setw(14) << d << endl;

  }
  
  mwrotation.close();
  
}
			 


int main () {
  
  // MW
  cout << "Initialising MW ... " << endl;
  string sph_cache_name_mw  =     "data/run068s22h/SLGridSph.cache.mw.run068s22h";
  string model_file_mw                      = "data/run068s22h/SLGridSph.mw.s22h";
  string coef_file_mw       = "data/run068s22h/simpleoutcoef.nofac.mw.run068s22h";
  string orient_file_mw             = "data/run068s22h/mw.orient.run068s22h.smth";

  SphExpansion* MW;
  MW = new SphExpansion(sph_cache_name_mw, model_file_mw, coef_file_mw, orient_file_mw);

  // LMC
  cout << "Initialising LMC ... " << endl;
  string sph_cache_name_lmc      = "data/run068s22h/SLGridSph.cache.lmc.run068s22h";
  string model_file_lmc                      = "data/run068s22h/SLGridSph.lmc.s22h";
  string coef_file_lmc       = "data/run068s22h/simpleoutcoef.nofac.lmc.run068s22h";
  string orient_file_lmc             = "data/run068s22h/lmc.orient.run068s22h.smth";

  SphExpansion* LMC;
  LMC = new SphExpansion(sph_cache_name_lmc, model_file_lmc, coef_file_lmc, orient_file_lmc);


  bool onlymonopole = true;
  array_type2 mwcoefs,lmccoefs;
  select_coefficient_time(0.0, MW->coeftable, mwcoefs);
  select_coefficient_time(0.0, LMC->coeftable, lmccoefs);

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
		      1000,
		      rotationfile, onlymonopole);

  
  rotationfile="tests/LMCrotation.txt";

  make_rotation_curve(LMC,
		      lmccoefs,
		      0.1,
		      120.,
		      1000,
		      rotationfile);
  

  vector<double> xinit(3);
  xinit[0] = 300.;
  xinit[1] = 0.;
  xinit[2] = 0.;

  vector<double> vinit(3);
  vinit[0] = 0.;
  vinit[1] = 220./1.4;
  //vinit[1] = 150./1.4;
  //vinit[1] = 0.;
  vinit[2] = 0.;

  array_type2 orbit;

  double dtintegrate = 0.01; // this is fine for circular orbits, needs to be ~0.003 for centre-crossing
  
  leapfrog(MW, mwcoefs,
	      xinit, vinit,
	      2000,
	      dtintegrate,
	      orbit);

  string orbitfile="tests/circularorbitMW.txt";
  print_orbit(orbit,orbitfile);

  // also print an LMC orbit
  xinit[0] = 100.;
  vinit[1] = 90.;

  leapfrog(LMC, lmccoefs,
	      xinit, vinit,
	      2000,
	      dtintegrate,
	      orbit);

  orbitfile="tests/circularorbitLMC.txt";
  print_orbit(orbit,orbitfile);
  
}
