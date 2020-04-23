/*
the main expansion example code, demonstrating the capability to find the forces at specific points.

compile string: 
clang++ -I/opt/local/include -L/opt/local/lib expansion.cc -o expansion

clean version, MSP 22 April 2020
revised for two-component models to spec, MSP 23 April 2020

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdio.h>

// boost includes
#include "boost/multi_array.hpp"

// MSP headers
// converters from r to xi (for mapping tables)
#include "scaling.h"

// some basic basis elements
#include "basis.h"

// basic transformations from cartesian to spherical (and back)
#include "transform.h"

// basic translations from virial to physical units (and back)
#include "translate.h"

// the cachefile stuff, also brings in the modelfile stuff
#include "sphcache.h"

// the coefficient stuff
#include "sphcoefs.h"

// the orientation stuff for centering the expansions
#include "sphorient.h"

using namespace std;

// create 2- and 3-d array types from boost
typedef boost::multi_array<double, 3> array_type3;
typedef boost::multi_array<double, 2> array_type2;

class SphExpansion
{
private:

  // doesn't need to be exposed to the outside world
  SphModel modeltable;
  
  void initialise(string sph_cache_name,
		  string model_file,
		  string coef_file,
		  string orient_file);

public:
  // the constructor
  SphExpansion(string sph_cache_name,
	       string model_file,
	       string coef_file,
	       string orient_file);

  // expose the important expansion data
  SphCache cachetable;
  SphOrient orient;
  SphCoefs coeftable;
  
  void get_pot_coefs(int l, int indx, int nmax, array_type2& coefs, array_type2& potd, array_type2& dpot, double *p, double *dp);

  void determine_fields_at_point_sph(SphCache& cachetable,
				     array_type2& coefs,
				     double r, double theta, double phi, 
				     double& potl0, double& potl, 
				     double& potr, double& pott, double& potp);

};

SphExpansion::SphExpansion(string sph_cache_name,
			   string model_file,
			   string coef_file,
			   string orient_file)
{
  initialise(sph_cache_name, model_file, coef_file, orient_file);
}

void SphExpansion::initialise(string sph_cache_name,
			      string model_file,
			      string coef_file,
			      string orient_file)
{
  // pull in the parts for the expansion
  read_sph_cache(sph_cache_name, SphExpansion::cachetable);

  read_model(model_file, SphExpansion::modeltable);

  read_coef_file (coef_file, SphExpansion::coeftable);

  read_orient (orient_file, SphExpansion::orient);
  
  // finish setting up the model
  init_table(SphExpansion::modeltable, SphExpansion::cachetable);

}



void SphExpansion::get_pot_coefs(int l, int indx, int nmax, array_type2& coefs, array_type2& potd, array_type2& dpot, double *p, double *dp)
{
  /*
    int l    : the harmonic order
    int indx : the indexed harmonic order (e.g. l and m values)
    
   */
  double pp, dpp;
  int i;

  pp = dpp = 0.0;

  for (i=0; i<nmax; i++) {
    pp  += potd[l][i] * coefs[indx][i];
    dpp += dpot[l][i] * coefs[indx][i];
  }

  *p = -pp;
  *dp = -dpp;
}



void SphExpansion::determine_fields_at_point_sph
(SphCache& cachetable,
 array_type2& coefs,
 double r, double theta, double phi, 
 double& potl0, double& potl, 
 double& potr, double& pott, double& potp)
{
  /*
  // skipping density for now --> decide later if interesting.
  
  see the equivalent exp call in SphericalBasis.cc

  */
  
  int l,loffset,moffset,m;
  double rs,fac1,fac2,fac3,fac4,costh,dp;
  double p,pc,dpc,ps,dps,dens;

  // block here, some problem with a zero in theta here. TBD.
  if (theta<1.e-6) theta = 1.e-6;
  costh = cos(theta);

  fac1 = 0.25/M_PI;

  array_type2 factrl;
  factorial(cachetable.LMAX, factrl);

  array_type2 legs, dlegs;
  dlegendre_R(cachetable.LMAX, costh, legs, dlegs);

  vector<double> cosm(cachetable.NMAX),sinm(cachetable.NMAX);
  sinecosine_R(cachetable.LMAX, phi, cosm, sinm);

  array_type2 potd,dpot;
  get_dpotl(r, cachetable, potd, dpot);

  // compute the monopole values
  get_pot_coefs(0, 0, cachetable.NMAX, coefs, potd, dpot, &p, &dp);
  potl = fac1*p;
  potr = fac1*dp;
  pott = potp = 0.0;

  // l loop
  for (l=1, loffset=1; l<=cachetable.LMAX; loffset+=(2*l+1), l++) {
    
    // m loop
    for (m=0, moffset=0; m<=l; m++) {
      fac1 = (2.0*l+1.0)/(4.0*M_PI);
      if (m==0) {
	fac2 = fac1*legs[l][m];
	
	get_pot_coefs(l, loffset+moffset, cachetable.NMAX, coefs, potd, dpot, &p, &dp);
	potl += fac2*p;
	potr += fac2*dp;
	pott += fac1*dlegs[l][m]*p;
	moffset++;
      }
      else {
	fac2 = 2.0 * fac1 * factrl[l][m];
	fac3 = fac2 *  legs[l][m];
	fac4 = fac2 * dlegs[l][m];
	
	get_pot_coefs(l, loffset+moffset,   cachetable.NMAX, coefs, potd, dpot, &pc, &dpc);
	get_pot_coefs(l, loffset+moffset+1, cachetable.NMAX, coefs, potd, dpot, &ps, &dps);
	
	potl += fac3*( pc*cosm[m] + ps*sinm[m]);
	potr += fac3*(dpc*cosm[m] + dps*sinm[m]);
	pott += fac4*( pc*cosm[m] +  ps*sinm[m]);
	potp += fac3*(-pc*sinm[m] +  ps*cosm[m])*m;
	
	moffset +=2;
      }
    }
  } 
}


// the best time is T=2.224: the present day. uniquely set per simulation
double reference_time = 2.224;


void return_forces_mw_and_lmc(SphExpansion* MW, SphExpansion* LMC,
			      array_type2 mwcoefs, array_type2 lmccoefs,
			      double t, double x, double y, double z,
			      double& fx, double& fy, double& fz)
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

  // reset time to have the correct system zero
  tvir += reference_time;

  // get the present-day MW coordinates: the zero of the system
  double zerox,zeroy,zeroz;
  zerox = MW->orient.xspline(reference_time);
  zeroy = MW->orient.yspline(reference_time);
  zeroz = MW->orient.zspline(reference_time);

  // grab the centres of the expansions at the specified times in exp reference space
  double xcenmw,ycenmw,zcenmw;
  xcenmw = MW->orient.xspline(tvir) - zerox;
  ycenmw = MW->orient.yspline(tvir) - zeroy;
  zcenmw = MW->orient.zspline(tvir) - zeroz;
  
  double xcenlmc,ycenlmc,zcenlmc;
  xcenlmc = LMC->orient.xspline(tvir) - zerox;
  ycenlmc = LMC->orient.yspline(tvir) - zeroy;
  zcenlmc = LMC->orient.zspline(tvir) - zeroz;

  // this would be a nice step to check the timing for. are we really hurting ourselves with the expensive spline interpolation?
  // i.e. should we set up a simple interpolator?
  //array_type2 coefsinmw,coefsinlmc;
  //select_coefficient_time(tvir, MW->coeftable, coefsinmw);
  //select_coefficient_time(tvir, LMC->coeftable, coefsinlmc);

  // the answer is yes, spline is expensive. moving the call outside helps; a more simple interpolation may help even further.
  // let's put a pin in this and if it's too slow, circle back.

  double rtmp,phitmp,thetatmp;
  double tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys;
  
  cartesian_to_spherical(xvir-xcenmw, yvir-ycenmw, zvir-zcenmw, rtmp, phitmp, thetatmp);
  
  MW->determine_fields_at_point_sph(MW->cachetable, mwcoefs,
				rtmp,thetatmp,phitmp,
				tpotl0,tpotl,
				fr,ft,fp);

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
				fr, fp, ft,
				fxtmp, fytmp, fztmp);

  //virial_to_physical_length(x,y,z,xphys,yphys,zphys);
  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  // same procedure for LMC
  cartesian_to_spherical(xvir-xcenlmc, yvir-ycenlmc, zvir-zcenlmc, rtmp, phitmp, thetatmp);
  
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
  cout << "The time is " << intime << " Gyr from pericentre." << endl;

  // get the time in exp system units by 1) converting the virial units, 2) applying time offset for present-day
  double tvir;
  physical_to_virial_time(intime,tvir);
  tvir += reference_time;

  // this would be a nice step to time. are we really hurting ourselves with the expensive spline interpolation?
  // i.e. should we set up a simple interpolator?
  array_type2 mwcoefs,lmccoefs;
  select_coefficient_time(tvir, MW->coeftable, mwcoefs);
  select_coefficient_time(tvir, LMC->coeftable, lmccoefs);
  // see discussion above in return_forces_mw_and_lmc
  
  double xin,yin,zin;
  double fx,fy,fz;

  // demo to make a rotation curve
  for (int xx=0; xx<100; xx++) {

    // the location in inertial space of the points to check (yin=zin=0, just checking x-axis right now)
    xin = xx*4. + 0.1; // in kpc

    return_forces_mw_and_lmc(MW, LMC,
			     mwcoefs, lmccoefs,
			     intime, xin, 0., 0.,
			     fx, fy, fz);

    cout << setw(14) << xin << setw(14) << fx << setw(14) << fy << setw(14) << fz << endl;

  }
  
}
