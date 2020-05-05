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
#include "cylcache.h"

//#include "rotate.h"

#include "cylcoefs.h"
#include "sphcoefs.h"

int main () {

  string cyl_cache_name_mw = "data/.eof.cache.001";
  string cyl_coef_name_mw = "data/outcoef.star.run001s";

  // pull in the parts for the MW disc cache
  //CylCache MWdisc;
  //read_cyl_cache(cyl_cache_name_mw, MWdisc);

  //array_type2 Vc,Vs;
  //get_pot(0.01, 0.0, MWdisc, Vc, Vs);

  //cout << Vc[0][0] << endl;

  CylCoefs MWdiscC;
  read_coef_file(cyl_coef_name_mw, MWdiscC);

  string sph_coef_name_mw = "data/outcoef.dark.run001s";
  SphCoefs MWhaloC;
  read_coef_file_raw(sph_coef_name_mw, MWhaloC);
  
}
