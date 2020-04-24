/*
tests for the initial MW and LMC
- generate a rotation curve
- integrate a circular orbit

compile string: 
clang++ -I/opt/local/include -L/opt/local/lib -Iinclude/ rotcurve.cc -o rotcurve

MSP 24 April 2020

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
			 string outfile)
{
  /*
  xmin and xmax come in PHYSICAL units
  and the rotation curve is returned in physical units
   */
  ofstream mwrotation;
  mwrotation.open(outfile);

  mwrotation << "# radius [kpc]; vcirc [km/s] ; f_x [km/s/s] ; f_y [km/s/s] ; f_z [km/s/s];" << endl;

  double dx = (xmax-xmin)/nsamples;
  double xin;
  double fx,fy,fz;

  // demo to make a rotation curve
  for (int xx=0; xx<nsamples; xx++) {

    // the location in inertial space of the points to check (yin=zin=0, just checking x-axis right now)
    xin = xx*dx + xmin; // in kpc

    return_forces(S,
		  coefs,
		  0.0, xin, 0., 0.,
		  fx, fy, fz);

    mwrotation << setw(14) << xin << setw(14) << sqrt(xin*fx) << setw(14) << fx << setw(14) << fy << setw(14) << fz << endl;

  }
  
  mwrotation.close();
  
}
			 


int main () {

  // obviously these would all be better read in . . . depends on how your code interfaces
  string sph_cache_name_mw = "data/SLGridSph.cache.mw.run068s10";
  string model_file_mw     = "data/SLGridSph.mw";
  string coef_file_mw      = "data/simpleoutcoef.nofac.mw.run068s10s";
  string orient_file_mw    = "data/mw.simpleorient.run068s10s";
 
  // pull in the parts for the MW
  SphExpansion* MW;
  MW = new SphExpansion(sph_cache_name_mw, model_file_mw, coef_file_mw, orient_file_mw);

  // try the LMC
  string sph_cache_name_lmc = "data/SLGridSph.cache.lmc.run068s10";
  string model_file_lmc     = "data/SLGridSph.lmc";
  string coef_file_lmc      = "data/simpleoutcoef.nofac.lmc.run068s10s";
  string orient_file_lmc    = "data/lmc.simpleorient.run068s10s";
 
  // pull in the parts for the LMC
  SphExpansion* LMC;
  LMC = new SphExpansion(sph_cache_name_lmc, model_file_lmc, coef_file_lmc, orient_file_lmc);

  
  array_type2 mwcoefs,lmccoefs;
  select_coefficient_time(0.0, MW->coeftable, mwcoefs);
  select_coefficient_time(0.0, LMC->coeftable, lmccoefs);
  
  
  string rotationfile="tests/MWrotation.txt";

  make_rotation_curve(MW,
		      mwcoefs,
		      0.1,
		      120.,
		      1000,
		      rotationfile);

  
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
  //vinit[1] = 220./1.4;
  //vinit[1] = 150./1.4;
  vinit[1] = 0.;
  vinit[2] = 0.;

  array_type2 orbit;
  
  leapfrog(MW, mwcoefs,
	      xinit, vinit,
	      1000,
	      0.003,
	      orbit);

  string orbitfile="tests/circularorbit.txt";

  print_orbit(orbit,orbitfile);
  
}
