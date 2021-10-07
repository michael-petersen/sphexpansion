/*
definitions for the SphExpansion class

MSP 24 Apr 2020 restructured
MSP  7 Oct 2020 fix density normalisation
MSP  4 Feb 2021 add dipole and quadrupole capability
MSP 25 May 2021 add ltrunc option; add self-gravitating coefficient calculation
MSP 28 Sep 2021 adjust radial order (nmax) truncation

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
  SphCoefs coeftable;
  
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
  SphCache cachetable; // does this actually have to be exposed?
  SphOrient orient;

  // expose the basic parameters of the expansion
  int LMAX; // maximum azimuthal order in expansion
  int NMAX; // maximum radial order in expansion
  int NUMT; // number of timesteps in the coefficients
  // can use these for checking that different files match!

  // get potential function weights
  void get_pot_coefs(int l, int indx, int nmax, array_type2& coefs, array_type2& potd, array_type2& dpot, double *p, double *dp);

  // get density function weights
  void get_dens_coefs(int l, int indx, int nmax, array_type2& coefs, array_type2& dend, double *dd);

  
  // the base spherical class
  // the flags at the end can specify which components you want to control specifically
  // monopole  =true means that only the monopole is considered, but at all spatial scales
  // dipole    =true means that only the monopole and the dipole are considered
  // quadrupole=true means that only the monopole and the quadrupole are considered
  // ltrunc          is an integer that specifies the maximum number of l harmonics to return
  //
  // notes: if monopole is true, dipole and quadrupole are forced to be false.
  //        you can set both dipole and quadrupole to true; you will get monopole+dipole+quadrupole.
  void determine_fields_at_point_sph(array_type2& coefs,
				     double r, double theta, double phi, 
				     double& potl0, double& potl, 
				     double& potr, double& pott,
				     double& potp,
				     bool monopole=false, bool dipole=false, bool quadrupole=false,
				     int ltrunc=1000);

  // version with density return
  void determine_fields_at_point_sph(array_type2& coefs,
				     double r, double theta, double phi, 
				     double& dens0, double& dens, 
				     double& potl0, double& potl, 
				     double& potr, double& pott,
				     double& potp,
				     bool monopole=false, bool dipole=false, bool quadrupole=false,
				     int ltrunc=1000);
  
  // cartesian forces wrapper function
  void return_forces(array_type2& coefs,
		     double x, double y, double z,
		     double& fx, double& fy, double& fz,
		     bool monopole=false, bool dipole=false, bool quadrupole=false,
		     int ltrunc=1000);

  // cartesian forces wrapper function
  void return_density(array_type2& coefs,
		      double x, double y, double z,
		      double& d,
		      bool monopole=false, bool dipole=false, bool quadrupole=false,
		      int ltrunc=1000);

  // coefficient interpolator
  //  this call may also be used as a coarse method to truncate coefficient series.
  //  however, it will not result in any speedup, as the coefficients are still evaluated (as zeros).
  void select_coefficient_time(double desired_time,
			       array_type2& coefs_at_time, int ntrunc=-1, int ltrunc=0);

  // return self-gravitating coefficients (i.e. apply spherical harmonic norm)
  void get_selfgravity_coefficients(array_type3& self_grav_coefs, bool monopolenorm=false);
  
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
  try {
    read_sph_cache(sph_cache_name, SphExpansion::cachetable);
  } catch (const char* msg) {
    cerr << msg << endl;
    exit(1);
  }

  read_model(model_file, SphExpansion::modeltable);

  read_coef_file (coef_file, SphExpansion::coeftable);
  LMAX = coeftable.LMAX;
  NMAX = coeftable.NMAX;
  NUMT = coeftable.NUMT;

  // if no orient file, assume zeros? 
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

void SphExpansion::get_dens_coefs(int l, int indx, int nmax, array_type2& coefs, array_type2& dend, double *dd)
{
  /*
    int l    : the harmonic order
    int indx : the indexed harmonic order (e.g. l and m values)
    
   */
  double daccum;
  int i;

  daccum = 0.0;

  for (i=0; i<nmax; i++)
    daccum  += dend[l][i] * coefs[indx][i];

  *dd = daccum;
}


void SphExpansion::determine_fields_at_point_sph
(array_type2& coefs,
 double r, double theta, double phi, 
 double& potl0, double& potl, 
 double& potr, double& pott, double& potp,
 bool monopole, bool dipole, bool quadrupole,
 int ltrunc)
{
  /*
  // version without density

  // no potl0 definition?
  
  see the equivalent exp call in SphericalBasis.cc

  */

  //cout << dipole << " " << quadrupole << endl;

  int numl = cachetable.LMAX;

  int l,loffset,moffset,m;
  double rs,fac1,fac2,fac3,fac4,costh,dp;
  double p,pc,dpc,ps,dps;//,dens;

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

    // advance loops if dipole or quadrupole flags are flown
    if (dipole && quadrupole && l>2) {
      continue;
    } else {
      if ( dipole && !quadrupole && l!=1) continue;
      if (!dipole &&  quadrupole && l!=2) continue;
    }

    // flag for selecting higher order terms
    if (l>ltrunc) continue;
    
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


void SphExpansion::determine_fields_at_point_sph(array_type2& coefs,
 double r, double theta, double phi, 
 double& dens0, double& dens, 
 double& potl0, double& potl, 
 double& potr, double& pott, double& potp,
 bool monopole, bool dipole, bool quadrupole,
 int ltrunc)
{
  /*
  // version WITH density
  
  // no potl0/dens0 definition?

  see the equivalent exp call in SphericalBasis.cc

  */


  int numl = cachetable.LMAX;
  
  int l,loffset,moffset,m;
  double rs,fac1,fac2,fac3,fac4,costh,dp;
  double p,pc,dpc,ps,dps,d,dc,ds;

  // block here, some problem with a zero in theta here. TBD.
  if (theta<1.e-6) theta = 1.e-6;
  costh = cos(theta);

  //fac0 = 4.*M_PI;
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
  dens = fac1*d;
  dens0 = fac1*d;


  // l loop
  if (monopole) return;

   array_type2 factrl;
  factorial(numl, factrl);

  array_type2 legs, dlegs;
  dlegendre_R(numl, costh, legs, dlegs);

  vector<double> cosm(cachetable.NMAX),sinm(cachetable.NMAX);
  sinecosine_R(numl, phi, cosm, sinm);
    
  for (l=1, loffset=1; l<=numl; loffset+=(2*l+1), l++) {

    // advance loops if dipole or quadrupole flags are flown
    if (dipole && quadrupole && l>2) {
      continue;
    } else {
      if ( dipole && !quadrupole && l!=1) continue;
      if (!dipole &&  quadrupole && l!=2) continue;
    }
    
    // flag for selecting higher order terms
    if (l>ltrunc) continue;
    
    // m loop
    for (m=0, moffset=0; m<=l; m++) {
      fac1 = (2.0*l+1.0)/(4.0*M_PI);
      if (m==0) {
	fac2 = fac1*legs[l][m];

	get_dens_coefs(l,loffset+moffset,cachetable.NMAX, coefs, dend, &d);
	dens += fac2*d;
	
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
	dens += fac3*(dc*cosm[m] + ds*sinm[m]);
	
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

void SphExpansion::return_forces(array_type2& coefs, double x, double y, double z,
				 double& fx, double& fy, double& fz,
				 bool monopole, bool dipole, bool quadrupole,
				 int ltrunc)
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
  
  determine_fields_at_point_sph(coefs, rtmp,thetatmp,phitmp,
				tpotl0,tpotl,
				fr,ft,fp,monopole,dipole,quadrupole,ltrunc);

  // DEEP debug
  //cout << setw(14) << rtmp << setw(14) << thetatmp << setw(14) << phitmp << setw(14) << fr << setw(14) << ft << setw(14) << fp << endl; 

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
				fr, fp, ft,
				fxtmp, fytmp, fztmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fx,fy,fz);

}


void SphExpansion::return_density(array_type2& coefs,
				  double x, double y, double z,
				  double& d,
				  bool monopole, bool dipole, bool quadrupole, int ltrunc)
{
  /*
    return density
   */

  // translate all times and positions into exp virial units
  double xvir,yvir,zvir;
  physical_to_virial_length(x,y,z,xvir,yvir,zvir);

  double rtmp,phitmp,thetatmp;
  double tpotl0,tpotl,fr,ft,fp,tdens0;
  d = 0;
  
  cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);
  
  determine_fields_at_point_sph(coefs,
				rtmp,thetatmp,phitmp,
				tdens0,d,
			        tpotl0,tpotl,
				fr,ft,fp,monopole,dipole,quadrupole,ltrunc);

#if DEEPDEBUGCOEFS
  cout << setw(14) << rtmp << setw(14) << thetatmp << setw(14) << phitmp << setw(14) << fr << setw(14) << ft << setw(14) << fp << endl; 
#endif

}

void SphExpansion::select_coefficient_time(double desired_time,
			     array_type2& coefs_at_time, int ntrunc, int ltrunc) {
  /*
    linear interpolation to get the coefficient matrix at a specific time

   time units must be virial time units to match the input coefficient
   table

   if ntrunc<0, will select all orders. if ntrunc>=0, will truncate at
   the specified n_max by setting all coefficients above order to 0.0.

   the ntrunc limit will be applied for all orders at or above ltrunc.
   by default, the ntruncation will be applied for all l orders.
   */

  int numl;

  numl = (coeftable.LMAX+1)*(coeftable.LMAX+1);

  coefs_at_time.resize(boost::extents[numl][coeftable.NMAX]);  
  
  // coeftable.t is assumed to be evenly spaced
  double dt = coeftable.t[1] - coeftable.t[0];

  int indx = (int)( (desired_time-coeftable.t[0])/dt);

  // guard against wanton extrapolation: should this stop the model?
  //if (indx<0) cerr << "select_coefficient_time: time prior to simulation start selected. setting to earliest step." << endl;

  // guard against going past the end of the simulation
  if (indx>coeftable.NUMT-2) cerr << "select_coefficient_time: time after to simulation end selected. setting to latest step." << endl;

#if DEEPDEBUGCOEFS
  cout << indx << endl;
#endif
  
  if (indx<0) {


    for (int l=0; l<numl; l++){
      for (int n=0; n<coeftable.NMAX; n++) {
        coefs_at_time[l][n] = coeftable.coefs[0][l][n];
      }
    }
    
  } else {

    // case where the simulations are not before the beginning of the simulation
    // interpolate from the two closest times

    double x1 = (coeftable.t[indx+1] - desired_time)/dt;
    double x2 = (desired_time - coeftable.t[indx])/dt;

#if DEEPDEBUGCOEFS
    cout << "dt=" << setw(16) << dt << "  t[indx+1]="  << setw(16) << coeftable.t[indx+1]
	                            << "  t[indx]="  << setw(16) << coeftable.t[indx]
                                    << "  desiredT="  << setw(16) << desired_time << endl;
    cout << "x1/x2=" << setw(16) << x1 << setw(14) << x2 << endl;
#endif

    for (int l=0; l<numl; l++){
      for (int n=0; n<coeftable.NMAX; n++) {
        coefs_at_time[l][n] = (x1*coeftable.coefs[indx][l][n] + x2*coeftable.coefs[indx+1][l][n]);
      }
    }
  }

  // go through and zero out non-selected orders
  //   note: will be applied for all l orders >= ltrunc

  if (ntrunc>=0) {
    for (int l=ltrunc*ltrunc; l<numl; l++){
    
      // start at order+1 so this call will work if n=0.
      for (int n=ntrunc+1; n<coeftable.NMAX; n++) {
        coefs_at_time[l][n] = 0.0;
      }
      
    }
  }
  
}



void SphExpansion::get_selfgravity_coefficients(array_type3& self_grav_coefs, bool monopolenorm)
{
  
  // function to produce coefficients that have been normalised with the spherical harmonic norm such that the self-gravity can be compared
  // this means including (1/(4*pi)) * (2*l+1) * ((l-m)!/(l+m)!)
  // this will set up another copy of the coefficients, so be warned: this could be large

  // The theory:
  // -The biorthogonality condition is the integral of the density and and the potential over 3d space.
  // -The inner product of \rho and \phi for the entire system is the *squared sum* of all the coefficients in the expansion.
  // -Physically, that is 2 times the gravitational potential energy.

  // Note that this function _does not_ return the squared sum, so if you want to compute energy, compare the squares of different terms.

  // todo: 
  // consider options for single lorder,morder return

  int numl = coeftable.LMAX;
  int numn = coeftable.NMAX;
  int l,loffset,moffset,m,n,t;

  self_grav_coefs.resize(boost::extents[coeftable.NUMT][(coeftable.LMAX+1)*(coeftable.LMAX+1)][coeftable.NMAX]);

  double fac1,fac2;
  double norm = 1.0;

  fac1 = 0.25/M_PI;
  
  array_type2 factrl;
  factorial(numl, factrl);

  for (t=0;t<coeftable.NUMT;t++) {

    // do the monopole
    for (n=0;n<numn;n++) self_grav_coefs[t][0][n] = fac1*coeftable.coefs[t][0][n];

    // if monopole norm, set up the norm
    if (monopolenorm) {
      norm = 0.0;

      // two choices of norms: the total monopole, or the lowest-order function only.
      
      // get the norm (the sum of the amplitude in all of the monopole terms, which hold the mass)
      //for (n=0;n<numn;n++) norm += fac1*coeftable.coefs[t][0][n];
      
      // get the simple norm (the power in the lowest-order function only)
      norm = fac1*coeftable.coefs[t][0][0];
      // now we are ready to norm the monopole
      for (n=0;n<numn;n++) self_grav_coefs[t][0][n] = self_grav_coefs[t][0][n]/norm;			     
    }

    // do the higher orders
    for (l=1, loffset=1; l<=numl; loffset+=(2*l+1), l++) {
      for (m=0, moffset=0; m<=l; m++) {
        fac1 = (2.0*l+1.0)/(4.0*M_PI);
        if (m==0) {
          for (n=0;n<numn;n++) self_grav_coefs[t][loffset+moffset][n] = (fac1/norm)*coeftable.coefs[t][loffset+moffset][n];
	  moffset++;
        } else {
	  fac2 = 2.0 * fac1 * factrl[l][m];
          for (n=0;n<numn;n++) self_grav_coefs[t][loffset+moffset  ][n] = (fac2/norm)*coeftable.coefs[t][loffset+moffset  ][n];
          for (n=0;n<numn;n++) self_grav_coefs[t][loffset+moffset+1][n] = (fac2/norm)*coeftable.coefs[t][loffset+moffset+1][n];
	  moffset += 2;
        }
      }  
    } // higher-order l loop

    
  } // time loop

}
