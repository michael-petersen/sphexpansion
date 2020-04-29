/*
example code:
integrate an orbit in the MW-LMC potential and compare to exp output

compile string: 
clang++ -I/opt/local/include -L/opt/local/lib -Iinclude/ fullintegrate.cc -o fullintegrate

MSP 28 Apr 2020 initial commit

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
#include "leapfrog.h"


// the best time is T=2.224: the present day. uniquely set per simulation.
double reference_time = 0;//2.224;




void return_forces_mw_and_lmc(SphExpansion* MW, SphExpansion* LMC,
			      array_type2 mwcoefs, array_type2 lmccoefs,
			      double t, double x, double y, double z,
			      double& fx, double& fy, double& fz, bool verbose)
{
  /*
    specs: take a time, x,y,z; return x,y,z forces, in physical units
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
    lmc_centre[j] -= zerocoords[j];
  }

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



void two_component_leapfrog(SphExpansion* MW,
			    array_type3 mwcoefs,
			    SphExpansion* LMC,
			    array_type3 lmccoefs,
			    vector<double> xinit,
			    vector<double> vinit,
			    int nint,
			    double dt,
			    array_type2& orbit)
{
  /*

   */
  double fx,fy,fz;

  
  double tvir;

  // include the forces for now
  orbit.resize(boost::extents[10][nint]);

  // initialise beginning values
  orbit[0][0] = xinit[0];
  orbit[1][0] = xinit[1];
  orbit[2][0] = xinit[2];
  orbit[3][0] = vinit[0];
  orbit[4][0] = vinit[1];
  orbit[5][0] = vinit[2];

  //now step forward one, using leapfrog (drift-kick-drift) integrator?
  //    https://en.wikipedia.org/wiki/Leapfrog_integration
  //
  int step = 1;

  // get the initial coefficient values: the time here is in tvir units, so always start with 0
  array_type2 tcoefsmw,tcoefslmc;
  select_coefficient_time(0., MW->coeftable, tcoefsmw);
  select_coefficient_time(0., LMC->coeftable, tcoefslmc);

  // not applying time offsets here; think about whether this is a problem
  double tphys;
  virial_to_physical_time(0.,tphys);
  // return forces for the initial step
  return_forces_mw_and_lmc(MW, LMC,
			     tcoefsmw, tcoefslmc,
			     tphys, orbit[0][0],orbit[1][0],orbit[2][0],
			     fx, fy, fz, false);

  orbit[6][0] = fx;
  orbit[7][0] = fy;
  orbit[8][0] = fz;

  int j;

  for (step=1; step<nint; step++) {

    // advance timestep: this is in physical units by definition.
    orbit[9][step] = dt*step;

    // find the current virial time
    physical_to_virial_time(dt*(step),tvir);

    // check that this is incrementing correctly? it is.
    //cout << "TVIR=" << tvir << endl;

    // get coefficients at the current virial time
    select_coefficient_time(tvir, MW->coeftable, tcoefsmw);
    select_coefficient_time(tvir, LMC->coeftable, tcoefslmc);

    // advance positions
    for (j=0; j<3; j++) {
      orbit[j][step] = orbit[j][step-1]   + (orbit[j+3][step-1]*dt  )  + (0.5*orbit[j+6][step-1]  * (dt*dt));
    }

    // calculate new forces: time goes in as physical time (e.g. kpc/km/s)
    return_forces_mw_and_lmc(MW, LMC,
			     tcoefsmw, tcoefslmc,
			     dt*(step-1), orbit[0][step],orbit[1][step],orbit[2][step],
			     fx, fy, fz, false);

    orbit[6][step] = fx;
    orbit[7][step] = fy;
    orbit[8][step] = fz;

    // advance velocities
    for (j=3; j<6; j++) {
      orbit[j][step] = orbit[j][step-1] + (0.5*(orbit[j+3][step-1]+orbit[j+3][step])  * dt );
    }
    
  }
  
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


  vector<double> zerocoords(3);
  return_centre(reference_time, MW->orient, zerocoords);

  cout << "Coordinate zero:" << setw(14) << zerocoords[0] << setw(14) << zerocoords[1] << setw(14) << zerocoords[2] << endl; 
  

  // get the time in exp system units by 1) converting the virial units, 2) applying time offset for present-day

  // initial position offset from centre
  // 0.25297400000000003 0.9793989999999999 -0.1661203
  vector<double> xphys(3),xinit(3);
  xphys[0] = 0.25297400000000003;
  xphys[1] = 0.9793989999999999;
  xphys[2] = -0.1661203;

  // or not transformed, 0.123143 0.558696 0.0459567
  xphys[0] = 0.123143;
  xphys[1] = 0.558696;
  xphys[2] = 0.0459567;

  // initial velocity as read from file
  // 0.457418 -0.464346 -0.124515
  vector<double> vxphys(3),vxinit(3);
  vxphys[0] = 0.457418;
  vxphys[1] = -0.464346;
  vxphys[2] = -0.124515;

  // 0.106459 0.0508656 0.103126
  xphys[0] = 0.106459;
  xphys[1] = 0.0508656;
  xphys[2] = 0.103126;
  //  -1.34426 1.69256 0.815283
  vxphys[0] = -1.34426;
  vxphys[1] = 1.69256;
  vxphys[2] = 0.815283;

  //0.567928 2.2214 -0.561674
  xphys[0] = 0.567928;
  xphys[1] = 2.2214;
  xphys[2] = -0.561674;
  //-0.00898586 -0.485074 -0.188528
  vxphys[0] = -0.00898586;
  vxphys[1] = -0.485074;
  vxphys[2] = -0.188528;

  
  virial_to_physical_length(xphys[0],xphys[1],xphys[2],xinit[0],xinit[1],xinit[2]);
  virial_to_physical_velocity(vxphys[0],vxphys[1],vxphys[2],vxinit[0],vxinit[1],vxinit[2]);

  cout << "Input pos/vel: " << xinit[0] << " " << xinit[1] << " " << xinit[2] << " " <<
    vxinit[0] << " " << vxinit[1] << " " << vxinit[2] << " " << endl;

  double nint=500;

  // call this time in kpc/km/s units
  double dt;
  // force sampling at the native rate of the exp simulation as an interpolation test
  virial_to_physical_time(0.001,dt);
  array_type2 orbit;

  two_component_leapfrog(MW, MW->coeftable.coefs, LMC, LMC->coeftable.coefs,xinit, vxinit, nint, dt, orbit);

  string orbitfile="tests/comparisonorbit6.txt";
  print_orbit(orbit,orbitfile);
  
}
