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
  //MW->get_selfgravity_coefficients(self_grav_coefs, true);
  //for (int i=0;i<MW->NUMT;i++) cout << self_grav_coefs[i][0][1] << endl;


  LMC->get_selfgravity_coefficients(self_grav_coefs, true, true);

  double tphys;
  double reference_time = 2.82;
  for (int i=0;i<LMC->NUMT;i++) {
    virial_to_physical_time(LMC->T[i], tphys);
    cout << LMC->T[i] << "-->" << tphys-reference_time << endl;
    //for (int j=0;j<(LMC->LMAX+1)*(LMC->LMAX+1);j++) cout << self_grav_coefs[i][j][0] << " ";
    for (int j=0;j<LMC->NMAX;j++) cout << self_grav_coefs[i][0][j] << " ";
    cout << endl;
  }

  cout << LMC->NUMT << endl;
  
  /*
  // example: keep 4 radial terms for l>1 (but keep all monopole terms)
  array_type2 trunc_coefs;
  MW->select_coefficient_time(0.0,trunc_coefs,4,1);
  for (int i=0;i<10;i++) {
    cout << "(0," << i << "):" << trunc_coefs[0][i] << endl;
    cout << "(1," << i << "):" << trunc_coefs[1][i] << endl;

  }
  */
}
