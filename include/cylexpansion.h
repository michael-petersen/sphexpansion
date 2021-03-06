/*
definitions for the CylExpansion class

MSP  5 May 2020 first commit
MSP 29 Sep 2020 test first model
MSP 13 Oct 2020 add monopole capability
MSP  4 Feb 2021 add dipole and quadrupole capability

this should be able to go much faster. obviously I'm passing too much
of something around, would like to track down what exactly is going on.

*/

#if STANDALONE
// converters from r to xi (for mapping tables)
#include "scaling.h"

// some basic basis elements
#include "basis.h"

// basic transformations from cartesian to spherical (and back)
#include "transform.h"

// basic translations from virial to physical units (and back)
#include "translate.h"

// the orientation stuff for centering the expansions
#include "sphorient.h"
#endif


// the cachefile stuff, also brings in the modelfile stuff
#include "cylcache.h"

// the coefficient stuff
#include "cylcoefs.h"




using namespace std;

// create 2- and 3-d array types from boost
typedef boost::multi_array<double, 4> array_type4;
typedef boost::multi_array<double, 3> array_type3;
typedef boost::multi_array<double, 2> array_type2;

class CylExpansion
{
private:

  CylCoefs coeftable;
  CylForce forcetable;
  CylCache cachetable;

  void initialise(string cyl_cache_name,
		  string coef_file,
		  string orient_file);
  
  
public:
  // the constructor
  CylExpansion(string cyl_cache_name,
	       string coef_file,
	       string orient_file);

  // expose the important expansion data: use the spherical code, it works here!
  SphOrient orient;

  // the base cylindrical class
  // the flags at the end can specify which components you want to control specifically
  // monopole  =true means that only the monopole is considered, but at all spatial scales
  // dipole    =true means that only the monopole and the dipole are considered
  // quadrupole=true means that only the monopole and the quadrupole are considered
  //
  // notes: if monopole is true, dipole and quadrupole are forced to be false.
  //        you can set both dipole and quadrupole to true; you will get monopole+dipole+quadrupole.
  void determine_fields_at_point_cyl(array_type2& coscoefs,
				     array_type2& sincoefs,
				     double r, double phi, double z, 
				     double& potl0, double& potl, 
				     double& potr, double& potp,
				     double& potz,
				     bool monopole=false, bool dipole=false, bool quadrupole=false);

  // cartesian forces wrapper function
  void return_forces(array_type2 coscoefs,
		     array_type2 sincoefs,
		     double x, double y, double z,
		     double& fx, double& fy, double& fz,
		     bool monopole=false, bool dipole=false, bool quadrupole=false);

  void select_coefficient_time(double desired_time,
			     array_type2& coscoefs_at_time,
			       array_type2& sincoefs_at_time);

  void get_table_forces(double r, double z, CylForce& forcetable);

};

CylExpansion::CylExpansion(string cyl_cache_name,
			   string coef_file,
			   string orient_file)
{
  initialise(cyl_cache_name, coef_file, orient_file);
}

void CylExpansion::initialise(string cyl_cache_name,
			      string coef_file,
			      string orient_file)
{
  // pull in the parts for the expansion
  read_cyl_cache(cyl_cache_name, CylExpansion::cachetable);

  read_coef_file (coef_file, CylExpansion::coeftable);

  read_orient (orient_file, CylExpansion::orient);

}




void CylExpansion::determine_fields_at_point_cyl(array_type2& coscoefs,
						 array_type2& sincoefs,
						 double r, double phi, double z, 
						 double& potl0, double& potl, 
						 double& fr, double& fp, double& fz,
						 bool monopole, bool dipole, bool quadrupole)
{
  /*
  // skipping density for now --> decide later if interesting.
  
  see the equivalent exp call in SphericalBasis.cc

  */
  
  int m,n;

  double ccos,ssin,fac;

  //get_table_forces(r, z, cachetable, forcetable);
  get_table_forces(r, z, forcetable);

  potl0  = 0.0;
  potl   = 0.0;
  fr     = 0.0;
  fz     = 0.0;
  fp     = 0.0;

  for (m=0; m<=cachetable.MMAX; m++) {

    if (monopole   &&   m>0) {
      return;
    }

    if (dipole && quadrupole) {
      if (m>2) continue;
    } else {
      
      if (dipole     && m!=0 && m!=1) {
        continue;
      }
    
      if (quadrupole && m!=0 && m!=2) {
        continue;
      }
    }

    ccos = cos(phi*m);
    ssin = sin(phi*m);

    for (n=0; n<cachetable.NORDER; n++) {

      fac = coscoefs[m][n] * ccos;

      if (m==0 && n==0) potl0 += fac * forcetable.potC[m][n];
      
      potl += fac * forcetable.potC[m][n];
      fr   += fac * forcetable.rforceC[m][n];
      fz   += fac * forcetable.zforceC[m][n];

      fac = coscoefs[m][n] * ssin;

      fp += fac * m * forcetable.potC[m][n];

      if (m) { // sine terms
	
	fac = sincoefs[m][n] * ssin;

	potl += fac * forcetable.potS[m][n];

	fr  += fac * forcetable.rforceS[m][n];
	fz  += fac * forcetable.zforceS[m][n];

	fac = -sincoefs[m][n] * ccos;

	fp  += fac * m * forcetable.potS[m][n];
	
      } // end sine loop
    } // end NORDER loop
  } // end MMAX loop

  
}

void CylExpansion::return_forces(array_type2 coscoefs,
		   array_type2 sincoefs,
		   double x, double y, double z,
		   double& fx, double& fy, double& fz,
		   bool monopole, bool dipole, bool quadrupole)
{
  /*
    force return from just one component, from the centre of the expansion
   */

  // translate all times and positions into exp virial units
  double xvir,yvir,zvir;
  physical_to_virial_length(x,y,z,xvir,yvir,zvir);

  double rtmp,phitmp;
  double potl0,potl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;
  
  cartesian_to_cylindrical(xvir, yvir, rtmp, phitmp);
  
  determine_fields_at_point_cyl(coscoefs, sincoefs,
				rtmp,phitmp,z,
				potl0,potl,
				fr,fp,fz,monopole,dipole,quadrupole);

  // DEEP debug
  //cout << setw(14) << rtmp << setw(14) << thetatmp << setw(14) << phitmp << setw(14) << fr << setw(14) << ft << setw(14) << fp << endl; 

  cylindrical_forces_to_cartesian(rtmp, phitmp,
				  fr, fp,
				  fxtmp, fytmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fx,fy,fz);

}



void CylExpansion::select_coefficient_time(double desired_time,
			     array_type2& coscoefs_at_time,
			     array_type2& sincoefs_at_time) {
  /*
    linear interpolation to get the coefficient matrix at a specific time

   time units must be virial time units to match the input coefficient table
   */

  // coeftable.t is assumed to be evenly spaced
  double dt = coeftable.t[1] - coeftable.t[0];

  int indx = (int)( (desired_time-coeftable.t[0])/dt);

  // guard against wanton extrapolation: should this stop the model?
  if (indx<0) cerr << "select_coefficient_time: time prior to simulation start selected. setting to earliest step." << endl;
  if (indx>coeftable.NUMT-2) cerr << "select_coefficient_time: time after to simulation end selected. setting to latest step." << endl;

  if (indx<0) indx = 0;
  if (indx>coeftable.NUMT-2) indx = coeftable.NUMT - 2;

  double x1 = (coeftable.t[indx+1] - desired_time)/dt;
  double x2 = (desired_time - coeftable.t[indx])/dt;

  coscoefs_at_time.resize(boost::extents[coeftable.MMAX+1][coeftable.NORDER]);
  sincoefs_at_time.resize(boost::extents[coeftable.MMAX+1][coeftable.NORDER]);

  for (int m=0; m<=coeftable.MMAX; m++){
    for (int n=0; n<coeftable.NORDER; n++) {
      coscoefs_at_time[m][n] = (x1*coeftable.coscoefs[indx][m][n] + x2*coeftable.coscoefs[indx+1][m][n]);

      if (m) sincoefs_at_time[m][n] = (x1*coeftable.sincoefs[indx][m][n] + x2*coeftable.sincoefs[indx+1][m][n]);
    }
  }
  
}


  
void CylExpansion::get_table_forces(double r, double z, CylForce& forcetable)
{

  // return 2d tables required to compute the forces
  
  forcetable.potC.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER]);
  forcetable.potS.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER]);

  forcetable.rforceC.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER]);
  forcetable.rforceS.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER]);

  forcetable.zforceC.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER]);
  forcetable.zforceS.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER]);


  if (z/cachetable.ASCALE > cachetable.Rtable) z =  cachetable.Rtable*cachetable.ASCALE;
  if (z/cachetable.ASCALE <-cachetable.Rtable) z = -cachetable.Rtable*cachetable.ASCALE;

  double X = (r_to_xi(r,cachetable.CMAP,cachetable.ASCALE) - cachetable.XMIN)/cachetable.dX;
  double Y = (z_to_y(z,cachetable.HSCALE) - cachetable.YMIN)/cachetable.dY;

  int ix = (int)X;
  int iy = (int)Y;

  if (ix < 0) {
    ix = 0;
  }
  if (iy < 0) {
    iy = 0;
  }
  
  if (ix >= cachetable.NUMX) {
    ix = cachetable.NUMX-1;
  }
  if (iy >= cachetable.NUMY) {
    iy = cachetable.NUMY-1;
  }

  double delx0 = (double)ix + 1.0 - X;
  double dely0 = (double)iy + 1.0 - Y;
  double delx1 = X - (double)ix;
  double dely1 = Y - (double)iy;

  double c00 = delx0*dely0;
  double c10 = delx1*dely0;
  double c01 = delx0*dely1;
  double c11 = delx1*dely1;

  for (int mm=0; mm<=cachetable.MMAX; mm++) {
    
    for (int n=0; n<cachetable.NORDER; n++) {

      forcetable.potC[mm][n] = 
	(
	 cachetable.potC[mm][n][ix  ][iy  ] * c00 +
	 cachetable.potC[mm][n][ix+1][iy  ] * c10 +
	 cachetable.potC[mm][n][ix  ][iy+1] * c01 +
	 cachetable.potC[mm][n][ix+1][iy+1] * c11 
	 );

      forcetable.rforceC[mm][n] = 
	(
	 cachetable.rforceC[mm][n][ix  ][iy  ] * c00 +
	 cachetable.rforceC[mm][n][ix+1][iy  ] * c10 +
	 cachetable.rforceC[mm][n][ix  ][iy+1] * c01 +
	 cachetable.rforceC[mm][n][ix+1][iy+1] * c11 
	 );

      forcetable.zforceC[mm][n] = 
	(
	 cachetable.zforceC[mm][n][ix  ][iy  ] * c00 +
	 cachetable.zforceC[mm][n][ix+1][iy  ] * c10 +
	 cachetable.zforceC[mm][n][ix  ][iy+1] * c01 +
	 cachetable.zforceC[mm][n][ix+1][iy+1] * c11 
	 );

      // get sine values for m>0
      if (mm) {

	forcetable.potS[mm][n] = 
	  (
	   cachetable.potS[mm][n][ix  ][iy  ] * c00 +
	   cachetable.potS[mm][n][ix+1][iy  ] * c10 +
	   cachetable.potS[mm][n][ix  ][iy+1] * c01 +
	   cachetable.potS[mm][n][ix+1][iy+1] * c11 
	   );

        forcetable.rforceS[mm][n] = 
	  (
	   cachetable.rforceS[mm][n][ix  ][iy  ] * c00 +
	   cachetable.rforceS[mm][n][ix+1][iy  ] * c10 +
	   cachetable.rforceS[mm][n][ix  ][iy+1] * c01 +
	   cachetable.rforceS[mm][n][ix+1][iy+1] * c11 
	   );

        forcetable.zforceS[mm][n] = 
	  (
	   cachetable.zforceS[mm][n][ix  ][iy  ] * c00 +
	   cachetable.zforceS[mm][n][ix+1][iy  ] * c10 +
	   cachetable.zforceS[mm][n][ix  ][iy+1] * c01 +
	   cachetable.zforceS[mm][n][ix+1][iy+1] * c11 
	   );
      }

    } // end NORDER loop
  } // end MMAX loop

}

