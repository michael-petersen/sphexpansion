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

 
  // open files to check trajectories
  string trajectoryfile="tests/trajectories.txt";
  ofstream traj;
  traj.open(trajectoryfile);

  traj << "# time ; mwx ; mwy ; mwz ; mwvx ; mwvy ; mwvz ; lmcx ; lmcy; lmcz; lmcvx; lmcvy; lmcvz;" << endl;

  // find the centres at the reference time
  vector<double> mw0centre(3),lmc0centre(3);
  vector<double> mw0velcentre(3),lmc0velcentre(3);

  double intime;
  
  for (int i=0; i<1000; i++) {

    intime = 0.004*i - 2.;
    return_centre(intime, LMC->orient, lmc0centre);
    return_centre(intime,  MW->orient,  mw0centre);
    return_vel_centre(intime, LMC->orient, lmc0velcentre);
    return_vel_centre(intime,  MW->orient,  mw0velcentre);

    traj << setw(14) << intime
	 << setw(14) << mw0centre[0]
	 << setw(14) << mw0centre[1]
	 << setw(14) << mw0centre[2]
	 << setw(14) << mw0velcentre[0]
	 << setw(14) << mw0velcentre[1]
	 << setw(14) << mw0velcentre[2]
	 << setw(14) << lmc0centre[0]
	 << setw(14) << lmc0centre[1]
	 << setw(14) << lmc0centre[2]
	 << setw(14) << lmc0velcentre[0]
	 << setw(14) << lmc0velcentre[1]
	 << setw(14) << lmc0velcentre[2]
	 << endl;
  


  }

      traj.close();
      
}

  
