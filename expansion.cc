/*
the main expansion example code, demonstrating the capability to find the forces at specific points.

compile string: 
clang++ -I/opt/local/include -L/opt/local/lib expansion.cc -o expansion

if interested, see the various testing codes in slreadtest.cc

clean version, MSP 22 April 2020

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

// basic transformations
#include "transform.h"

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
typedef boost::multi_array<double, 1> array_type1;

typedef vector<double> dbl_vector;
typedef boost::multi_array<vector<dbl_vector>, 2> efarray;


class SphExpansion
{
private:
  int LMAX;
  
  void initialise(string sph_cache_name,
			   string model_file,
			   string coef_file,
		  string orient_file,
		  SphCache& cachetable,
		  SphModel& modeltable,
		  SphCoefs& coeftable,
		  SphOrient& orient);
public:
  // the constructor
  SphExpansion(string sph_cache_name,
	       string model_file,
	       string coef_file,
	       string orient_file,
	       SphCache& cachetable,
			   SphModel& modeltable,
			   SphCoefs& coeftable,
			   SphOrient& orient);

  void get_pot_coefs(int l, int indx, int nmax, array_type2& coefs, array_type2& potd, array_type2& dpot, double *p, double *dp);

  void determine_fields_at_point_sph
(SphCache& cachetable,
 array_type2& coefs,
 double r, double theta, double phi, 
 double& potl0, double& potl, 
 double& potr, double& pott, double& potp);

  
  
};


SphExpansion::SphExpansion(string sph_cache_name,
			   string model_file,
			   string coef_file,
			   string orient_file,
			   SphCache& cachetable,
			   SphModel& modeltable,
			   SphCoefs& coeftable,
			   SphOrient& orient)
{
  initialise(sph_cache_name, model_file, coef_file, orient_file, cachetable, modeltable, coeftable,orient);
}

void SphExpansion::initialise(string sph_cache_name,
			   string model_file,
			   string coef_file,
			      string orient_file,
			      SphCache& cachetable,
			      SphModel& modeltable,
			      SphCoefs& coeftable,
			      SphOrient& orient)
{
  // pull in the parts for the MW
  //SphCache cachetable;
  read_sph_cache(sph_cache_name, cachetable);

  //SphModel modeltable;
  read_model(model_file, modeltable);

  //SphCoefs coeftable; 
  read_coef_file (coef_file, coeftable);

  // finish setting up the model
  init_table(modeltable, cachetable);

  //SphOrient orient;
  read_orient (orient_file, orient);
  
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


int main () {

  // obviously these would all be better read in...depends on how your code interfaces
  string sph_cache_name_mw = "data/SLGridSph.cache.mw.run068s10";
  string model_file_mw     = "data/SLGridSph.mw";
  string coef_file_mw      = "data/simpleoutcoef.nofac.mw.run068s10";
  string orient_file_mw    = "data/mw.simpleorient.run068s10";
 
  // pull in the parts for the MW
  SphCache cachetablemw;
  SphModel modeltablemw;
  SphCoefs coeftablemw; 
  SphOrient orientmw;
  SphExpansion* MW;
  MW = new SphExpansion(sph_cache_name_mw, model_file_mw, coef_file_mw, orient_file_mw, cachetablemw, modeltablemw, coeftablemw, orientmw);

  // pick a timestep to analyse from the coefs
  int timestep = 0;

  double intime = coeftablemw.t[timestep];
  cout << "The time is " << intime << endl;
    
  array_type2 coefsin = coeftablemw.coefs[timestep];

  // get the offsets from the inertial
  double xcen,ycen,zcen;
  xcen = orientmw.xspline(intime);
  ycen = orientmw.yspline(intime);
  zcen = orientmw.zspline(intime);

  // print out the inertial offsets
  cout << setw(18) << xcen << setw(18) << ycen << setw(18) << zcen << endl; 

  double tpotl0,tpotl,fr,ft,fp,fx,fy,fz;
  double xin,yin,zin,rin,phiin,thetain;
  
  // demo to make a rotation curve
  for (int xx=0; xx<200; xx++) {

  yin = xx*0.02 + 0.0001;
  
  cartesian_to_spherical(xin, yin, zin, rin, phiin, thetain);
  
  MW->determine_fields_at_point_sph(cachetablemw, coefsin,
				rin,thetain,phiin,
				tpotl0,tpotl,
				fr,ft,fp);

  spherical_forces_to_cartesian(rin, phiin, thetain,
				   fr, fp, ft,
				fx, fy, fz);

  cout << setw(18) << rin << setw(18) << fx << setw(18) << fy << setw(18) << fz << endl;

  }
  
}
