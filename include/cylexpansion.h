/*
definitions for the CylExpansion class

MSP 5 May 2020 first commit

this should be able to go much faster. obviously I'm passing too much
of something around, would like to track down what exactly is going on.

*/



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
#include "cylcache.h"

// the coefficient stuff
#include "cylcoefs.h"

// the orientation stuff for centering the expansions
#include "sphorient.h"

using namespace std;

// create 2- and 3-d array types from boost
typedef boost::multi_array<double, 4> array_type4;
typedef boost::multi_array<double, 3> array_type3;
typedef boost::multi_array<double, 2> array_type2;

class CylExpansion
{
private:
  
  void initialise(string cyl_cache_name,
		  string coef_file,
		  string orient_file);
  
  // read up on 'protected'
  
public:
  // the constructor
  CylExpansion(string cyl_cache_name,
	       string coef_file,
	       string orient_file);

  // expose the important expansion data
  CylCache cachetable;
  SphOrient orient;
  CylCoefs coeftable;

  // the base spherical class
  void determine_fields_at_point_cyl(CylCache& cachetable,
				     array_type2& coscoefs,
				     array_type2& sincoefs,
				     double r, double phi, double z, 
				     double& potl0, double& potl, 
				     double& potr, double& potp,
				     double& potz);

  // cartesian forces wrapper function
  void return_forces(CylExpansion* C,
		     array_type2 coscoefs,
		     array_type2 sincoefs,
		     double x, double y, double z,
		     double& fx, double& fy, double& fz);

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

  //read_orient (orient_file, CylExpansion::orient);

}




void CylExpansion::determine_fields_at_point_cyl
        (CylCache& cachetable,
         array_type2& coscoefs,
         array_type2& sincoefs,
         double r, double phi, double z, 
         double& potl0, double& potl, 
         double& fr, double& fp, double& fz)
{
  /*
  // skipping density for now --> decide later if interesting.
  
  see the equivalent exp call in SphericalBasis.cc

  */
  
  int m,n;

  double ccos,ssin,fac;
  
  CylForce forcetable;
  get_table_forces(r, z, cachetable, forcetable);

  potl0  = 0.0;
  potl   = 0.0;
  fr     = 0.0;
  fz     = 0.0;
  fp     = 0.0;

  for (m=0; m<=cachetable.MMAX; m++) {

    ccos = cos(phi*m);
    ssin = sin(phi*m);

    for (n=0; n<cachetable.NORDER; n++) {

      fac = coscoefs[m][n] * ccos;

      potl += fac * forcetable.potC[m][n];

      if (m==0 && n==0) potl0 += fac * forcetable.potC[m][n];

      fr   += fac * forcetable.rforceC[m][n];
	
      fz   += fac * forcetable.zforceC[m][n];

      //d += np.sum(fac * (densC[mm,:,ix,iy]*c00 + densC[mm,:,ix+1,iy  ]*c10 + densC[mm,:,ix,iy+1]*c01 + densC[mm,:,ix+1,iy+1]*c11));

      fac = coscoefs[m][n] * ssin;

      fp += fac * m * forcetable.potC[m][n];

      if (m) { // sine terms
	
	fac = sincoefs[m][n] * ssin;

	potl += fac * forcetable.potS[m][n];

	fr  += fac * forcetable.rforceS[m][n];

	fz  += fac * forcetable.zforceS[m][n];

	//d += np.sum(fac * ( densS[mm,:,ix,iy] * c00 + densS[mm,:,ix+1,iy  ] * c10 + densS[mm,:,ix,iy+1] * c01 + densS[mm,:,ix+1,iy+1] * c11 ));

	fac = -sincoefs[m][n] * ccos;

	fp  += fac * m * forcetable.potS[m][n];
	
      } // end sine loop
    } // end NORDER loop
  } // end MMAX loop

  
}

void CylExpansion::return_forces(CylExpansion* C,
		   array_type2 coscoefs,
		   array_type2 sincoefs,
		   double x, double y, double z,
		   double& fx, double& fy, double& fz)
{
  /*
    test force return from just one component, from the centre of the expansion
   */

  // translate all times and positions into exp virial units
  double xvir,yvir,zvir;
  physical_to_virial_length(x,y,z,xvir,yvir,zvir);

  double rtmp,phitmp;
  double potl0,potl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;
  
  cartesian_to_cylindrical(xvir, yvir, rtmp, phitmp);
  
  C->determine_fields_at_point_cyl(C->cachetable, coscoefs, sincoefs,
				   rtmp,phitmp,z,
				   potl0,potl,
				   fr,fp,fz);

  // DEEP debug
  //cout << setw(14) << rtmp << setw(14) << thetatmp << setw(14) << phitmp << setw(14) << fr << setw(14) << ft << setw(14) << fp << endl; 

  cylindrical_forces_to_cartesian(rtmp, phitmp,
				  fr, fp,
				  fxtmp, fytmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fx,fy,fz);

}
