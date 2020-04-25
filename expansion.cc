/*
the main expansion example code:
compute the forces at a given cartesian point resulting from the MW and the LMC

compile string: 
clang++ -I/opt/local/include -L/opt/local/lib -Iinclude/ expansion.cc -o expansion

MSP 22 Apr 2020 clean version
MSP 23 Apr 2020 revised for two-component models to spec

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


// the best time is T=2.224: the present day. uniquely set per simulation
double reference_time = 2.224;


void return_forces_mw_and_lmc(SphExpansion* MW, SphExpansion* LMC,
			      array_type2 mwcoefs, array_type2 lmccoefs,
			      double t, double x, double y, double z,
			      double& fx, double& fy, double& fz, bool verbose)
{
  /*
    specs: take a time, x,y,z; return x,y,z forces.
    input/output units must be physical
    where x0,y0,z0 = present-day galactic centre
       and      t0 = present day (so previous times are negative)

   if we don't want to pass the entire SphExpansion objects, the necessary pieces can be broken out in a fairly straightforward way.

   */
  
  // zero out forces
  fx = fy = fz = 0.;

  // translate all times and positions into exp virial units
  double tvir,xvir,yvir,zvir;
  physical_to_virial_time(t,tvir);
  physical_to_virial_length(x,y,z, xvir,yvir,zvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  //cout << "TIME=" << tvir << endl;

  // initialise the centre vectors
  vector<double> zerocoords(3),mw_centre(3),lmc_centre(3);

  // get the present-day MW coordinates: the zero of the system
  return_centre(reference_time, MW->orient, zerocoords);

  //cout << "Coordinate zero:" << setw(14) << zerocoords[0] << setw(14) << zerocoords[1] << setw(14) << zerocoords[2] << endl; 

  // get the centres of the expansions at the specified times in exp reference space
  return_centre(tvir,  MW->orient,  mw_centre);
  return_centre(tvir, LMC->orient, lmc_centre);

  // shift the expansion centres to the pericentre coordinate system
  for (int j=0;j<=2;j++) {
    mw_centre[j]  -= zerocoords[j];
    //cout << setw(14) << mw_centre[j];
    lmc_centre[j] -= zerocoords[j];
    //cout << setw(14) << lmc_centre[j];
  }
  //cout << endl;

  if (verbose) {
    cout << "MW virial centre (x,y,z)=(" << mw_centre[0] << ","<< mw_centre[1] << ","<< mw_centre[2] << ")" <<endl;
    cout << "LMC virial centre (x,y,z)=(" << lmc_centre[0] << ","<< lmc_centre[1] << ","<< lmc_centre[2] << ")" <<endl;
  }

  double rtmp,phitmp,thetatmp;
  double tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys;

  // compute spherical coordinates in the frame of the MW expansion
  cartesian_to_spherical(xvir-mw_centre[0], yvir-mw_centre[1], zvir-mw_centre[2], rtmp, phitmp, thetatmp);

  //cout << setw(14) << rtmp << setw(14) << phitmp << setw(14) << thetatmp << endl;
  
  MW->determine_fields_at_point_sph(MW->cachetable, mwcoefs,
				rtmp,thetatmp,phitmp,
				tpotl0,tpotl,
				fr,ft,fp);

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
				fr, fp, ft,
				fxtmp, fytmp, fztmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  // same procedure for LMC
  cartesian_to_spherical(xvir-lmc_centre[0], yvir-lmc_centre[1], zvir-lmc_centre[2], rtmp, phitmp, thetatmp);

  //cout << setw(14) << rtmp << setw(14) << phitmp << setw(14) << thetatmp << endl;

  LMC->determine_fields_at_point_sph(LMC->cachetable, lmccoefs,
				rtmp,thetatmp,phitmp,
				tpotl0,tpotl,
				fr,ft,fp);

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
				fr, fp, ft,
				fxtmp, fytmp, fztmp);

  // reset to physical units
  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add LMC force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;
  
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

  // set a specific timestep here.
  double intime = -4.225; // in Gyr. -4.225 is the earliest recorded time in the model

  if (intime<0)  cout << "The time is " << -intime << " Gyr before pericentre." << endl;
  if (intime>=0)  cout << "The time is " << intime << " Gyr after pericentre." << endl;

  // get the time in exp system units by 1) converting the virial units, 2) applying time offset for present-day
  double tvir;
  physical_to_virial_time(intime,tvir);
  tvir += reference_time;

  array_type2 mwcoefs,lmccoefs;
  select_coefficient_time(tvir, MW->coeftable, mwcoefs);
  select_coefficient_time(tvir, LMC->coeftable, lmccoefs);
  
  double xin,yin,zin;
  double fx,fy,fz;

  ofstream combinedforces;
  combinedforces.open("tests/forcecurvex.txt");

  combinedforces << "# x [kpc] ; f_x [km/s/s] ; f_y [km/s/s] ; f_z [km/s/s];" << endl;
  // demo to make a force curve
  for (int xx=0; xx<100; xx++) {

    // the location in inertial space of the points to check (yin=zin=0, checking x-axis right now)
    xin = xx*4. + 0.1; // in kpc

    return_forces_mw_and_lmc(MW, LMC,
			     mwcoefs, lmccoefs,
			     intime, xin, 0., 0.,
			     fx, fy, fz, false);

    
    combinedforces << setw(14) << xin << setw(14) << fx << setw(14) << fy << setw(14) << fz << endl;

  }

  combinedforces.close();
  
}
