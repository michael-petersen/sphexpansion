/*
example code:
extract self-gravitating coefficients

compile string: 
clang++ -I/opt/local/include -L/opt/local/lib -I../include/ extractcoefs.cc -o obj/extractcoefs

MSP 25 May 2021 basic implementation

*/

#define STANDALONE 1


#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdio.h>

// boost includes
#include "boost/multi_array.hpp"

// the spherical expansion headers
#include "expansion.h"

// the cylindrical expansion headers
#include "cylexpansion.h"



int main () {

  // MW
  cout << "Initialising MW ... " << endl;
  string sph_cache_name_mw  = "/Volumes/External1/Disk076/SLGridSph.mw.run7mld";
  string model_file_mw      = "/Volumes/External1/Disk076/ErkalMW.model";
  string coef_file_mw       = "/Volumes/External1/Disk076/simpleoutcoef.nofac.mw.run7mld";
  string orient_file_mw     = "/Volumes/External1/Disk076/mw.orient.run7mld.smth";

  SphExpansion* MW;
  MW = new SphExpansion(sph_cache_name_mw, model_file_mw, coef_file_mw, orient_file_mw);

  // LMC
  cout << "Initialising LMC ... " << endl;
  string sph_cache_name_lmc = "/Volumes/External1/Disk076/SLGridSph.lmc.run7mld";
  string model_file_lmc     = "/Volumes/External1/Disk076/ErkalLMC.model";
  string coef_file_lmc      = "/Volumes/External1/Disk076/simpleoutcoef.nofac.lmc.run7mld";
  string orient_file_lmc    = "/Volumes/External1/Disk076/lmc.orient.run7mld.smth";

  SphExpansion* LMC;
  LMC = new SphExpansion(sph_cache_name_lmc, model_file_lmc, coef_file_lmc, orient_file_lmc);

  array_type3 self_grav_coefs;
  MW->get_selfgravity_coefficients(self_grav_coefs, -1,-1, true);
  for (int i=1;i<1000;i++)
  cout << self_grav_coefs[i][16][0] << endl;
  
}
