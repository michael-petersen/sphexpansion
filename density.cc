/*
tests for the initial MW and LMC
- density under deforming conditions
- made to track down a bug.

compile string: 
clang++ -I/opt/local/include -L/opt/local/lib -Iinclude/ density.cc -o density

MSP 13 Oct 2020 first written.

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

void make_density(SphExpansion* S1,
		  array_type2 coefs1,
		  SphExpansion* S2,
		  array_type2 coefs2,
		  double reftime,
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

  vector<double> mwcentre(3),lmccentre(3);
  return_centre(reftime,  S1->orient,  mwcentre);
  double mwxvp,mwyvp,mwzvp;
  virial_to_physical_length(mwcentre[0],mwcentre[1],mwcentre[2],mwxvp,mwyvp,mwzvp);

  return_centre(reftime,  S2->orient, lmccentre);
  double lmcxvp,lmcyvp,lmczvp;
  virial_to_physical_length(lmccentre[0],lmccentre[1],lmccentre[2],lmcxvp,lmcyvp,lmczvp);

  cout << setw(14) <<  mwxvp << setw(14) <<  mwyvp << setw(14) <<  mwzvp << endl;
  cout << setw(14) << lmcxvp << setw(14) << lmcyvp << setw(14) << lmczvp << endl;

  if (monopole) cout << "Only using monopole for rotation curve." << endl;
  
  ofstream densityfile;
  densityfile.open(outfile);

  densityfile << "#  y [kpc] ; z [km/s/s] ; mwdensity [Msun/pc^3] ; lmcdensity [Msun/pc^3];" << endl;

  double dx = (xmax-xmin)/nsamples;
  double zin,yin;
  double mwd,lmcd,mwphysdens,lmcphysdens;

  // demo to make a rotation curve
  for (int xx=0; xx<nsamples; xx++) {

    // the location in inertial space of the points to check (xin=0, just checking y-z plane right now)
    zin = xx*dx + xmin; // in kpc

    for (int yy=0; yy<nsamples; yy++) {

      yin = yy*dx + xmin;

      
      S1->return_density(S1,
			 coefs1,
			 0., yin, zin,
		         mwd, false);
      
      S2->return_density(S2,
			 coefs2,
			 0.-lmcxvp, yin-lmcyvp, zin-lmczvp,
		         lmcd, false);
    
      virial_to_physical_density( mwd,  mwphysdens);
      virial_to_physical_density(lmcd, lmcphysdens);
    
      densityfile  << setw(14) << yin        << setw(14) << zin <<
                      setw(14) << mwphysdens << setw(14) << lmcphysdens << endl;

    }
  }
  
  densityfile.close();
  
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

  double reftime = 1.50;
  bool onlymonopole = false;

  array_type2 mwcoefs,lmccoefs;
  select_coefficient_time(reftime, MW->coeftable, mwcoefs);
  select_coefficient_time(reftime, LMC->coeftable, lmccoefs);
  
  string dfile="tests/MW_LMC_density.txt";
  make_density(MW,mwcoefs,LMC,lmccoefs,reftime,-100,100.,100,dfile, false);


  
}
