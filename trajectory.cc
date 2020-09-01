/*
tests to compare the trajectories between two different components

compile string: 
clang++ -I/opt/local/include -L/opt/local/lib -Iinclude/ trajectory.cc -o trajectory

MSP 18 Aug 2020 first version

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

  /*
  m2 is 20%: reference_time   = 1.47787
  m2p5 is 25%: reference_time = 1.45625s
  m3 is 30%: reference_time   = 1.4355
   */
  double reference_time = 1.50;

  // find the centres at the reference time
  vector<double> mw0centre(3),lmc0centre(3);
  return_centre(reference_time, LMC->orient, lmc0centre);
  return_centre(reference_time,  MW->orient,  mw0centre);

  cout << "MW C.O.M. x=" << mw0centre[0] << " y=" <<
    mw0centre[1] << " z=" << mw0centre[2] << endl;
  
  cout << "LMC C.O.M. x=" << lmc0centre[0] << " y=" <<
    lmc0centre[1] << " z=" << lmc0centre[2] << endl;

  // find the velocity centres at the reference time
  vector<double> mw0velcentre(3),lmc0velcentre(3),lmcvelcentre(3);
  return_vel_centre(reference_time, LMC->orient, lmc0velcentre);
  return_vel_centre(reference_time,  MW->orient,  mw0velcentre);

  cout << "MW C.O.V. vx=" << mw0velcentre[0] << " vy=" <<
    mw0velcentre[1] << " vz=" << mw0velcentre[2] << endl;
  
  cout << "LMC C.O.V. vx=" << lmc0velcentre[0] << " vy=" <<
    lmc0velcentre[1] << " vz=" << lmc0velcentre[2] << endl;

    return_vel_centre(0.0, MW->orient, lmcvelcentre);
  cout << "MW T=0. C.O.V. vx=" << lmcvelcentre[0] << " vy=" <<
    lmcvelcentre[1] << " vz=" << lmcvelcentre[2] << endl;
}
