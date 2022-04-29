/*
definitions for the CylExpansion class

MSP  5 May 2020 first commit
MSP 29 Sep 2020 test first model
MSP 13 Oct 2020 add monopole capability
MSP  4 Feb 2021 add dipole and quadrupole capability
MSP 22 Dec 2021 fix spacing queries in select_coefficient_time

this should be able to go much faster. obviously I'm passing too much
of something around, would like to track down what exactly is going on.

*/
#ifndef CYLEXPANSION_H
#define CYLEXPANSION_H

// important preprocessor flags
#include "flags.h"

// converters from r to xi (for mapping tables)
#include "common/scaling.h"

// some basic basis elements
#include "common/basis.h"

// basic transformations from cartesian to spherical (and back)
#include "common/transform.h"

// basic translations from virial to physical units (and back)
#include "common/translate.h"

// the orientation stuff for centering the expansions
#include "common/orient.h"

// the cachefile stuff, also brings in the modelfile stuff
#include "cylinder/cylcache.h"

// the coefficient stuff
#include "cylinder/cylcoefs.h"

int CYLHARMONICDEFAULT = 2047;

//using namespace std;
using std::cout, std::cerr, std::endl, std::setw, std::vector, std::ifstream, std::ios, std::string, std::ofstream, std::istringstream;

// Eigen MatrixXd, std::vector <MatrixXd>
#include <Eigen/StdVector>
#include <Eigen/Dense>
using Eigen::MatrixXd;


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
  // behaviour for flags is now overridden by harmonicflag, which uses basis.check_flags as a binary bit flag to determine which orders to set. 2047 (default) will enable all values at l<=10
  void determine_fields_at_point_cyl(MatrixXd& coscoefs,
                                     MatrixXd& sincoefs,
                                     double r, double phi, double z,
                                     double& potl0, double& potl,
                                     double& potr, double& potp,
                                     double& potz, int harmonicflag=CYLHARMONICDEFAULT);

  // cartesian forces wrapper function
  void return_forces(MatrixXd coscoefs,
                     MatrixXd sincoefs,
                     double x, double y, double z,
                     double& fx, double& fy, double& fz, int harmonicflag=CYLHARMONICDEFAULT);

  void select_coefficient_time(double desired_time,
                               MatrixXd& coscoefs_at_time,
                               MatrixXd& sincoefs_at_time);

  void get_table_forces(double r, double z, CylForce& forcetable);

}; // end class definition

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
  try {
    read_cyl_cache(cyl_cache_name, CylExpansion::cachetable);
  } catch (const char* msg) {
    cerr << msg << endl;
    exit(1);
  }

  try {
    read_coef_file(coef_file, CylExpansion::coeftable);
  } catch (const char* msg) {
    cerr << msg << endl;
    exit(1);
  }

  try {
    read_orient   (orient_file, CylExpansion::orient);
  } catch (const char* msg) {
    cerr << msg << endl;
    exit(1);
  }

}




void CylExpansion::determine_fields_at_point_cyl(MatrixXd& coscoefs,
                                                 MatrixXd& sincoefs,
                                                 double r, double phi, double z,
                                                 double& potl0, double& potl,
                                                 double& fr, double& fp, double& fz,
                                                 int harmonicflag)
{
  /*
  @IMPROVE: no density call available here.

  */

  double ccos,ssin,fac;

  get_table_forces(r, z, forcetable);

  potl0  = 0.0;
  potl   = 0.0;
  fr     = 0.0;
  fz     = 0.0;
  fp     = 0.0;

  for (int m=0; m<=cachetable.MMAX; m++) {

    // check harmonic flag before proceeding
    if (check_flags(harmonicflag,m)==0) continue;

    ccos = cos(phi*m);
    ssin = sin(phi*m);

    for (int n=0; n<cachetable.NORDER; n++) {

      fac = coscoefs(m,n) * ccos;

      if (m==0 && n==0) potl0 += fac * forcetable.potC(m,n);

      potl += fac * forcetable.potC(m,n);
      fr   += fac * forcetable.rforceC(m,n);
      fz   += fac * forcetable.zforceC(m,n);

      fac = coscoefs(m,n) * ssin;

      fp += fac * m * forcetable.potC(m,n);

      if (m) { // sine terms

        fac = sincoefs(m,n) * ssin;

        potl += fac * forcetable.potS(m,n);

        fr  += fac * forcetable.rforceS(m,n);
        fz  += fac * forcetable.zforceS(m,n);

        fac = -sincoefs(m,n) * ccos;

        fp  += fac * m * forcetable.potS(m,n);

      } // end sine loop
    } // end NORDER loop
  } // end MMAX loop


}

void CylExpansion::return_forces(MatrixXd coscoefs,
                                 MatrixXd sincoefs,
                                 double x, double y, double z,
                                 double& fx, double& fy, double& fz,
                                 int harmonicflag)
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
  fr,fp,fz,harmonicflag);

#if DEEPDEBUGCOEFS
  std::cout << "cylexpansion.h::return_forces:" << setw(14) << rtmp << setw(14) << phitmp << setw(14) << z << setw(14) << fr << setw(14) << fp << setw(14) << fz << std::endl;
#endif

  cylindrical_forces_to_cartesian(rtmp, phitmp,
                                  fr, fp,
                                  fxtmp, fytmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fx,fy,fz);

}



void CylExpansion::select_coefficient_time(double desired_time,
                                           MatrixXd& coscoefs_at_time,
                                           MatrixXd& sincoefs_at_time)
{
  /*
    linear interpolation to get the coefficient matrix at a specific time

   time units must be virial time units to match the input coefficient table

   @IMPROVE: could add harmonic flag here to speed up calculation if not using higher orders
   */

  // starting at the first indx, stop when we get to the matching time
  int indx = 0;
  while (coeftable.t[indx]<=desired_time) {
    indx ++;
  }

  // reset by one
  indx --;

  // guard against wanton extrapolation: should this stop the model?
  // indx<0 is not possible in this construction
  // if (indx<0) cerr << "select_coefficient_time: time prior to simulation start selected. setting to earliest step." << endl;
  if (indx>coeftable.NUMT-2) std::cerr << "select_coefficient_time: time after to simulation end selected. setting to latest step." << std::endl;

  if (indx<0) indx = 0;
  if (indx>coeftable.NUMT-2) indx = coeftable.NUMT - 2;

  // check the local spacing on coeftable.t (can be globally nonuniform)
  double dt = coeftable.t[indx+1] - coeftable.t[indx];

#if DEEPDEBUGTIME
  std::cout << "cylexpansion.h::select_coefficient_time:" << setw(14) << desired_time << setw(14) << coeftable.t[0] << setw(14) << coeftable.t[indx] << setw(14) << dt << std::endl;
#endif

  double x1 = (coeftable.t[indx+1] - desired_time)/dt;
  double x2 = (desired_time - coeftable.t[indx])/dt;

  // does this resizing hurt us? we should only ever actually have to do it once.
  coscoefs_at_time.resize(coeftable.MMAX+1,coeftable.NORDER);
  sincoefs_at_time.resize(coeftable.MMAX+1,coeftable.NORDER);

  coscoefs_at_time = (x1*coeftable.coscoefs[indx] + x2*coeftable.coscoefs[indx+1]);
  sincoefs_at_time = (x1*coeftable.sincoefs[indx] + x2*coeftable.sincoefs[indx+1]);
  //for (int m=0; m<=coeftable.MMAX; m++){
  //  for (int n=0; n<coeftable.NORDER; n++) {
  //    coscoefs_at_time(m,n) = (x1*coeftable.coscoefs[indx](m,n) + x2*coeftable.coscoefs[indx+1](m,n));

  //    if (m) sincoefs_at_time(m,n) = (x1*coeftable.sincoefs[indx](m,n) + x2*coeftable.sincoefs[indx+1](m,n));
  //  }
  //}

}



void CylExpansion::get_table_forces(double r, double z, CylForce& forcetable)
{

  // return 2d tables required to compute the forces

  forcetable.potC.resize(coeftable.MMAX+1,coeftable.NORDER);
  forcetable.potS.resize(coeftable.MMAX+1,coeftable.NORDER);

  forcetable.rforceC.resize(coeftable.MMAX+1,coeftable.NORDER);
  forcetable.rforceS.resize(coeftable.MMAX+1,coeftable.NORDER);

  forcetable.zforceC.resize(coeftable.MMAX+1,coeftable.NORDER);
  forcetable.zforceS.resize(coeftable.MMAX+1,coeftable.NORDER);


  if (z/cachetable.ASCALE > cachetable.Rtable) z =  cachetable.Rtable*cachetable.ASCALE;
  if (z/cachetable.ASCALE <-cachetable.Rtable) z = -cachetable.Rtable*cachetable.ASCALE;

  double X = (r_to_xi_cyl(r,cachetable.CMAPR,cachetable.ASCALE) - cachetable.XMIN)/cachetable.dX;
  double Y = (z_to_y_cyl(z,cachetable.CMAPZ,cachetable.HSCALE) - cachetable.YMIN)/cachetable.dY;

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

      forcetable.potC(mm,n) =
        (
         cachetable.potC[mm][n](ix  ,iy  ) * c00 +
         cachetable.potC[mm][n](ix+1,iy  ) * c10 +
         cachetable.potC[mm][n](ix  ,iy+1) * c01 +
         cachetable.potC[mm][n](ix+1,iy+1) * c11
         );

      forcetable.rforceC(mm,n) =
        (
         cachetable.rforceC[mm][n](ix  ,iy  ) * c00 +
         cachetable.rforceC[mm][n](ix+1,iy  ) * c10 +
         cachetable.rforceC[mm][n](ix  ,iy+1) * c01 +
         cachetable.rforceC[mm][n](ix+1,iy+1) * c11
         );

      forcetable.zforceC(mm,n) =
        (
         cachetable.zforceC[mm][n](ix  ,iy  ) * c00 +
         cachetable.zforceC[mm][n](ix+1,iy  ) * c10 +
         cachetable.zforceC[mm][n](ix  ,iy+1) * c01 +
         cachetable.zforceC[mm][n](ix+1,iy+1) * c11
         );

      // get sine values for m>0
      if (mm) {

        forcetable.potS(mm,n) =
          (
           cachetable.potS[mm][n](ix  ,iy  ) * c00 +
           cachetable.potS[mm][n](ix+1,iy  ) * c10 +
           cachetable.potS[mm][n](ix  ,iy+1) * c01 +
           cachetable.potS[mm][n](ix+1,iy+1) * c11
           );

        forcetable.rforceS(mm,n) =
          (
           cachetable.rforceS[mm][n](ix  ,iy  ) * c00 +
           cachetable.rforceS[mm][n](ix+1,iy  ) * c10 +
           cachetable.rforceS[mm][n](ix  ,iy+1) * c01 +
           cachetable.rforceS[mm][n](ix+1,iy+1) * c11
           );

        forcetable.zforceS(mm,n) =
          (
           cachetable.zforceS[mm][n](ix  ,iy  ) * c00 +
           cachetable.zforceS[mm][n](ix+1,iy  ) * c10 +
           cachetable.zforceS[mm][n](ix  ,iy+1) * c01 +
           cachetable.zforceS[mm][n](ix+1,iy+1) * c11
           );
      }

    } // end NORDER loop
  } // end MMAX loop

}

#endif
