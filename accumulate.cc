/*
example code:
accumulate the coefficients from a particle distribution

compile string: 
clang++ -I/opt/local/include -L/opt/local/lib -Iinclude/ accumulate.cc -o accumulate

MSP 30 Apr 2020 initial commit

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdio.h>

// boost includes
#include "boost/multi_array.hpp"

// the expansion headers
#include "expansion.h"

// for print_orbit only, other integration pieces are below
#include "accumulate.h"



int main () {

  // obviously these would all be better read in . . . depends on how your code interfaces
  string sph_cache_name_mw = "data/SLGridSph.cache.mw.run068s10";
  string model_file_mw     = "data/SLGridSph.mw";
  string coef_file_mw      = "data/simpleoutcoef.nofac.mw.run068s10s";
  string orient_file_mw    = "data/mw.simpleorient.run068s10s";
 
  // pull in the parts for the MW
  SphExpansion* MW;
  MW = new SphExpansion(sph_cache_name_mw, model_file_mw, coef_file_mw, orient_file_mw);

}
