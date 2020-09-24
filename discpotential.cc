/*
example code:
compute forces from a disc

compile string: 
clang++ -I/opt/local/include -L/opt/local/lib -Iinclude/ discpotential.cc -o discpotential

MSP 5 May 2020 initial commit

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdio.h>

// boost includes
#include "boost/multi_array.hpp"

// the expansion headers
//#include "expansion.h"

// for print_orbit only, other integration pieces are below
//#include "cylcache.h"

//#include "rotate.h"

//#include "cylcoefs.h"
//#include "sphcoefs.h"

// the class for computing cylindrical forces
#include "cylexpansion.h"

//#include <yaml-cpp/yaml.h>	      // YAML support

void make_rotation_curve(CylExpansion* C,
			 array_type2 coscoefs,
			 array_type2 sincoefs,
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
  mwrotation.open(outfile, ios::out);

  mwrotation << "# radius [kpc]; vcirc [km/s] ; f_x [km/s/s] ; f_y [km/s/s] ; f_z [km/s/s];" << endl;

  double dx = (xmax-xmin)/nsamples;
  double xin;
  double fx,fy,fz;

  // demo to make a rotation curve
  for (int xx=0; xx<nsamples; xx++) {


    // the location in inertial space of the points to check (yin=zin=0, just checking x-axis right now)
    xin = xx*dx + xmin; // in kpc

    cout << xin << endl;

    C->return_forces(C,
		  coscoefs,
		  sincoefs,
		  xin, 0., 0.,
		  fx, fy, fz);

    mwrotation << setw(14) << xin << setw(14) << sqrt(xin*-fx) << setw(14) << fx << setw(14) << fy << setw(14) << fz << endl;

  }
  
  mwrotation.close();
  
}
			 



int main () {

  string cyl_cache_name_mw = "data/.eof.cache.001";
  string cyl_coef_name_mw = "data/outcoef.star.run001s";
  //string cyl_orient_name_mw = "data/outcoef.star.run001s";
  string cyl_orient_name_mw = "data/run068s22h/mw.orient.run068s22h.smth";

  CylExpansion* MWDisc;
  MWDisc = new CylExpansion(cyl_cache_name_mw, cyl_coef_name_mw, cyl_orient_name_mw);

  cout << "Step 1" << endl;

  array_type2 mwcoscoefs,mwsincoefs;
  select_coefficient_time(0.0, MWDisc->coeftable, mwcoscoefs, mwsincoefs);

  cout << "Step 2" << endl;
    
  string rotationfile="tests/MWdiscrotation.txt";

  make_rotation_curve(MWDisc,
		      mwcoscoefs,
		      mwsincoefs,
		      0.1,
		      30.,
		      100,
		      rotationfile);
  
}
