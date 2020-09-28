/*
definitions for the SphExpansion class

MSP 24 Apr 2020 restructured

*/

// turn off the inclusion of boilerplate stuff for CylExpansion
#undef STANDALONE
#define STANDALONE 0

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

  // get potential function weights
  void get_pot_coefs(int l, int indx, int nmax, array_type2& coefs, array_type2& potd, array_type2& dpot, double *p, double *dp);

  // get density function weights
  void get_dens_coefs(int l, int indx, int nmax, array_type2& coefs, array_type2& dend, double *p);

  
  // the base spherical class
  void determine_fields_at_point_sph(SphCache& cachetable,
				     array_type2& coefs,
				     double r, double theta, double phi, 
				     double& potl0, double& potl, 
				     double& potr, double& pott,
				     double& potp, bool monopole=false);

  // version with density return
  void determine_fields_at_point_sph(SphCache& cachetable,
				     array_type2& coefs,
				     double r, double theta, double phi, 
				     double& dens0, double& dens, 
				     double& potl0, double& potl, 
				     double& potr, double& pott,
				     double& potp, bool monopole=false);
  
  // cartesian forces wrapper function
  void return_forces(SphExpansion* S,
		     array_type2 coefs,
		     double x, double y, double z,
		     double& fx, double& fy, double& fz,
		     bool monopole=false);

  // cartesian forces wrapper function
  void return_density(SphExpansion* S,
		     array_type2 coefs,
		     double x, double y, double z,
		     double& d,
		     bool monopole=false);

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

void SphExpansion::get_dens_coefs(int l, int indx, int nmax, array_type2& coefs, array_type2& dend, double *d)
{
  /*
    int l    : the harmonic order
    int indx : the indexed harmonic order (e.g. l and m values)
    
   */
  double dd;
  int i;

  dd = 0.0;

  for (i=0; i<nmax; i++)
    dd  += dend[l][i] * coefs[indx][i];

  *d = dd;
}


void SphExpansion::determine_fields_at_point_sph
(SphCache& cachetable,
 array_type2& coefs,
 double r, double theta, double phi, 
 double& potl0, double& potl, 
 double& potr, double& pott, double& potp, bool monopole)
{
  /*
  // version without density

  // no potl0 definition?
  
  see the equivalent exp call in SphericalBasis.cc

  */

  int numl;
  if (monopole) {
    numl = 1;
  } else {
    numl = cachetable.LMAX;
  }
  
  int l,loffset,moffset,m;
  double rs,fac1,fac2,fac3,fac4,costh,dp;
  double p,pc,dpc,ps,dps,dens;

  // block here, some problem with a zero in theta here. TBD.
  if (theta<1.e-6) theta = 1.e-6;
  costh = cos(theta);

  fac1 = 0.25/M_PI;

  array_type2 potd,dpot;
  get_dpotl(r, cachetable, potd, dpot);

  // is this ever evaluating the l=0,n>0 terms??

  // compute the monopole values
  get_pot_coefs(0, 0, cachetable.NMAX, coefs, potd, dpot, &p, &dp);
  potl = fac1*p;
  potr = fac1*dp;
  pott = potp = 0.0;

  // l loop
  if (monopole) return;

   array_type2 factrl;
  factorial(numl, factrl);

  array_type2 legs, dlegs;
  dlegendre_R(numl, costh, legs, dlegs);

  vector<double> cosm(cachetable.NMAX),sinm(cachetable.NMAX);
  sinecosine_R(numl, phi, cosm, sinm);
    
  for (l=1, loffset=1; l<=numl; loffset+=(2*l+1), l++) {
    
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

  // zero out pott/potp if below resolution limit for forces
  // (e.g. centre crossing problems)
  if (r<cachetable.RMIN) {
    pott = 0.;
    potp = 0.;
  }
  
  
}


void SphExpansion::determine_fields_at_point_sph
(SphCache& cachetable,
 array_type2& coefs,
 double r, double theta, double phi, 
 double& dens0, double& dens, 
 double& potl0, double& potl, 
 double& potr, double& pott, double& potp, bool monopole)
{
  /*
  // version WITH density
  
  // no potl0/dens0 definition?

  see the equivalent exp call in SphericalBasis.cc

  */

  int numl;
  if (monopole) {
    numl = 1;
  } else {
    numl = cachetable.LMAX;
  }
  
  int l,loffset,moffset,m;
  double rs,fac1,fac2,fac3,fac4,costh,dp;
  double p,pc,dpc,ps,dps,d,dc,ds;

  // block here, some problem with a zero in theta here. TBD.
  if (theta<1.e-6) theta = 1.e-6;
  costh = cos(theta);

  fac1 = 0.25/M_PI;

  array_type2 potd,dpot,dend;
  get_dpotl_density(r, cachetable, potd, dpot, dend);

  // is this ever evaluating the l=0,n>0 terms??

  // compute the monopole values
  get_pot_coefs(0, 0, cachetable.NMAX, coefs, potd, dpot, &p, &dp);
  potl = fac1*p;
  potr = fac1*dp;
  pott = potp = 0.0;

  get_dens_coefs(0, 0, cachetable.NMAX, coefs, dend, &d);
  dens = d*fac1*fac1;
  dens0 = d*fac1*fac1;



  // l loop
  if (monopole) return;

   array_type2 factrl;
  factorial(numl, factrl);

  array_type2 legs, dlegs;
  dlegendre_R(numl, costh, legs, dlegs);

  vector<double> cosm(cachetable.NMAX),sinm(cachetable.NMAX);
  sinecosine_R(numl, phi, cosm, sinm);
    
  for (l=1, loffset=1; l<=numl; loffset+=(2*l+1), l++) {
    
    // m loop
    for (m=0, moffset=0; m<=l; m++) {
      fac1 = (2.0*l+1.0)/(4.0*M_PI);
      if (m==0) {
	fac2 = fac1*legs[l][m];
	
	get_dens_coefs(l,loffset+moffset,cachetable.NMAX, coefs, dend ,&d);
	dens += fac1*fac2*d;
	
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

	get_dens_coefs(l,loffset+moffset,  cachetable.NMAX, coefs, dend, &dc);
	get_dens_coefs(l,loffset+moffset+1,cachetable.NMAX, coefs, dend, &ds);
	dens += fac1*fac3*(pc*cosm[m] + ps*sinm[m]);
	
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

  // zero out pott/potp if below resolution limit for forces
  // (e.g. centre crossing problems)
  if (r<cachetable.RMIN) {
    pott = 0.;
    potp = 0.;
    dens = cachetable.d0[0]; // set to the smallest value of density. check if this is a good idea
  }
  
  
}

void SphExpansion::return_forces(SphExpansion* S,
		   array_type2 coefs,
		   double x, double y, double z,
		   double& fx, double& fy, double& fz,
		   bool monopole)
{
  /*
    test force return from just one component, from the centre of the expansion
   */

  // translate all times and positions into exp virial units
  double xvir,yvir,zvir;
  physical_to_virial_length(x,y,z,xvir,yvir,zvir);

  double rtmp,phitmp,thetatmp;
  double tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;
  
  cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);
  
  S->determine_fields_at_point_sph(S->cachetable, coefs,
				rtmp,thetatmp,phitmp,
				tpotl0,tpotl,
				   fr,ft,fp,monopole);

  // DEEP debug
  //cout << setw(14) << rtmp << setw(14) << thetatmp << setw(14) << phitmp << setw(14) << fr << setw(14) << ft << setw(14) << fp << endl; 

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
				fr, fp, ft,
				fxtmp, fytmp, fztmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fx,fy,fz);

}


void SphExpansion::return_density(SphExpansion* S,
		   array_type2 coefs,
		   double x, double y, double z,
		   double& d,
		   bool monopole)
{
  /*
    return density
   */

  // translate all times and positions into exp virial units
  double xvir,yvir,zvir;
  physical_to_virial_length(x,y,z,xvir,yvir,zvir);

  double rtmp,phitmp,thetatmp;
  double tpotl0,tpotl,fr,ft,fp,tdens0;
  
  cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);
  
  S->determine_fields_at_point_sph(S->cachetable, coefs,
				   rtmp,thetatmp,phitmp,
				   tdens0,d,
				   tpotl0,tpotl,
				   fr,ft,fp,monopole);

  // DEEP debug
  //cout << setw(14) << rtmp << setw(14) << thetatmp << setw(14) << phitmp << setw(14) << fr << setw(14) << ft << setw(14) << fp << endl; 


}
