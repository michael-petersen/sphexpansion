/*
header file for fullintegrate.cc capabilities

MSP 21 Apr 2022 cleaned version v0.2.3

*/
#ifndef FULLINTEGRATE_H
#define FULLINTEGRATE_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdio.h>

// eigen includes
#include <Eigen/Dense>
using Eigen::MatrixXd;

// set model parameters
#include "modelfiles.h"

// the spherical expansion headers
#include "sphexpansion.h"

// the cylindrical expansion headers
#include "cylexpansion.h"


class MWLMC
{
private:

  void initialise();

public:

  // the expansion classes: these might be better private
  SphExpansion* MW;
  SphExpansion* LMC;
  CylExpansion* MWD;

  // the constructor (no type, no arguments, all defined in modelfiles.h)
  MWLMC();

  // return all fields for the MW halo (in the frame of the MW halo)
  std::vector<double> mwhalo_fields(double t, double x, double y, double z,
                                    bool globalframe=false,
                                    int mwhharmonicflag=127, bool verbose=false);

  // return all fields for the LMC halo
  std::vector<double> lmc_fields(double t, double x, double y, double z,
                                 bool globalframe=false,
                                 int lmcarmonicflag=127, bool verbose=false);

  // return all fields for the MW disc
  // NOTE: density does not work here. not enabled yet for cylindrical expansions.
  //       leave as an inspirational placeholder
  std::vector<double>  mwd_fields(double t, double x, double y, double z,
                                  bool globalframe=false,
                                  int mwdharmonicflag=127, bool verbose=false);

  // return total forces
  // @IMPROVE: write all potential as well
  std::vector<double> all_forces(double t, double x, double y, double z,
                                 bool globalframe=true,
                                 int mwhharmonicflag=127, int mwdharmonicflag=127, int lmcharmonicflag=127,
                                 bool verbose=false);

  // version for integration, to avoid extra allocations
  void all_forces_coefs(MatrixXd mwcoefs, MatrixXd lmccoefs, MatrixXd mwdcoscoefs, MatrixXd mwdsincoefs,
                        double t, double x, double y, double z,
                        double& fx, double& fy, double& fz,
                        int mwhharmonicflag=127, int mwdharmonicflag=127, int lmcharmonicflag=127,
                        bool verbose=false);

  // if integrating just the MW
  void mw_forces_coefs(MatrixXd mwcoefs,
                       MatrixXd mwdcoscoefs, MatrixXd mwdsincoefs,
                       double t, double x, double y, double z,
                       double& fx, double& fy, double& fz,
                       int mwhharmonicflag, int mwdharmonicflag,
                       bool verbose=false);

  // compute an orbit integration in all three components
  // this uses PHYSICAL units by nature
  // lengths in kpc
  // velocities in km/s
  // time in Gyr
  MatrixXd orbit(vector<double> xinit,
                 vector<double> vinit,
                 double tbegin=-2.5,
                 double tend=0.0,
                 double dt=0.002,
                 int mwhharmonicflag=127, int mwdharmonicflag=127, int lmcharmonicflag=127,
                 bool discframe=true);

  // compute an orbit integration using only the initial Milky Way potential
  MatrixXd mworbit(vector<double> xinit,
                   vector<double> vinit,
                   double tbegin=-2.5,
                   double tend=0.0,
                   double dt=0.002);

   // compute an orbit rewind in all three components
   MatrixXd rewind(vector<double> xinit,
                   vector<double> vinit,
                   double dt=0.002,
                   int mwhharmonicflag=127, int mwdharmonicflag=127, int lmcharmonicflag=127,
                   double rewind_time=2.5,
                   bool discframe = true);

  // get centres of the expansions with PHYSICAL time input
  std::vector<double> get_expansion_centres_physical(double t, bool verbose=false);

  // get centres of the expansions with VIRIAL time input
  std::vector<double> get_expansion_centres_virial(double tvir, bool verbose=false);

  std::vector<double> get_lmc_centre_virial(double tvir, bool verbose);

  MatrixXd get_trajectories(double dt=native_timestep, bool virial=false);

  MatrixXd get_lmc_trajectory(double dt=native_timestep);

  // print an orbit array
  void print_orbit(MatrixXd orbit, string orbitfile);
};


MWLMC::MWLMC()
{
  initialise();
};

void MWLMC::initialise()
{
  // MW
  cout << "Initialising MW ... " << endl;
  MW = new SphExpansion(sph_cache_name_mw, model_file_mw, coef_file_mw, orient_file_mw);

  // LMC
  cout << "Initialising LMC ... " << endl;
  LMC = new SphExpansion(sph_cache_name_lmc, model_file_lmc, coef_file_lmc, orient_file_lmc);

  // MW
  cout << "Initialising MW disc ... " << endl;
  MWD = new CylExpansion(cyl_cache_name_mw, cyl_coef_name_mw, cyl_orient_name_mw);
};


void MWLMC::print_orbit(MatrixXd orbit, string orbitfile)
{
  ofstream outorbit;
  outorbit.open(orbitfile);

  outorbit << "# t [Gyr]; x [kpc]; y [kpc]; z [kpc]; vx [km/s] ; vy [km/s] ; vz [km/s] ; f_x [km/s/s] ; f_y [km/s/s] ; f_z [km/s/s];" << endl;

  for (int i=0; i<orbit.cols(); i++) {

    outorbit << setw(14) << orbit(9,i);

    for (int j=0; j<9; j++) {
      outorbit << setw(14) << orbit(j,i);
    }
    outorbit << endl;
  }

  outorbit.close();

}

std::vector<double> MWLMC::mwhalo_fields(double t, double x, double y, double z,
                                         bool globalframe,
                                         int mwhharmonicflag, bool verbose)
{
  /*
    always comes out in the frame of the expansion: translate to inertial before input, if desired

   */

   double tvir,xvir,yvir,zvir,xtmp,ytmp,ztmp;
   vector<double> zerocoords(3);

   if (globalframe) {

     // get the present-day MWD coordinates: the zero of the system
     return_centre(reference_time, MWD->orient, zerocoords);

     // shift the expansion centres to the pericentre coordinate system
     xtmp = x - zerocoords[0];
     ytmp = y - zerocoords[1];
     ztmp = z - zerocoords[2];

   } else {
     xtmp = x;
     ytmp = y;
     ztmp = z;
   }

  // translate all times and positions into exp virial units
  physical_to_virial_time(t,tvir);
  physical_to_virial_length(xtmp,ytmp,ztmp,xvir,yvir,zvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  // get coefficients
  MatrixXd mwcoefs;
  MW->select_coefficient_time(tvir, mwcoefs);

  double rtmp,phitmp,thetatmp,r2tmp;
  double dens0,denstmp,tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys,pphys,dphys;

  // compute spherical coordinates in the frame of the MW expansion
  cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);

  // get all field values
  MW->determine_fields_at_point_sph(mwcoefs,
                                    rtmp,thetatmp,phitmp,
                                    dens0,denstmp,
                                    tpotl0,tpotl,
                                    fr,ft,fp,
                                    mwhharmonicflag);

  // convert to cartesian
  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  // translate all quantities to physical units
  virial_to_physical_density(denstmp,dphys);
  virial_to_physical_potential(tpotl,pphys);
  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // initialise output array
  std::vector<double> output(5);

  output[0] = fxphys;
  output[1] = fyphys;
  output[2] = fzphys;
  output[3] = dphys;
  output[4] = pphys;

  return output;

}


std::vector<double> MWLMC::lmc_fields(double t, double x, double y, double z,
                                      bool globalframe,
                                      int lmcharmonicflag, bool verbose)
{
  /*
    always comes out in the frame of the expansion: translate to inertial before input, if desired

   */
   // translate all times and positions into exp virial units

  // translate all times and positions into exp virial units
  double tvir,xvir,yvir,zvir;
  physical_to_virial_time(t,tvir);
  physical_to_virial_length(x,y,z, xvir,yvir,zvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  MatrixXd lmccoefs;
  LMC->select_coefficient_time(tvir, lmccoefs);

  double rtmp,phitmp,thetatmp,r2tmp;
  double dens0,denstmp,tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys,pphys,dphys;

  // compute spherical coordinates in the frame of the LMC expansion
  cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);

  // get all field values
  LMC->determine_fields_at_point_sph(lmccoefs,
                                     rtmp,thetatmp,phitmp,
                                     dens0,denstmp,
                                     tpotl0,tpotl,
                                     fr,ft,fp,
                                     lmcharmonicflag);

  // convert to cartesian
  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  // translate all quantities to physical units
  virial_to_physical_density(denstmp,dphys);
  virial_to_physical_potential(tpotl,pphys);
  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // initialise output array
  std::vector<double> output(5);

  output[0] = fxphys;
  output[1] = fyphys;
  output[2] = fzphys;
  output[3] = dphys;
  output[4] = pphys;

  return output;

}


std::vector<double> MWLMC::mwd_fields(double t, double x, double y, double z,
                                      bool globalframe,
                                      int mwdharmonicflag,
                                      bool verbose)
{
  /*

   */

  // translate all times and positions into exp virial units
  double tvir,xvir,yvir,zvir;
  physical_to_virial_time(t,tvir);
  physical_to_virial_length(x,y,z, xvir,yvir,zvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  // get coefficients
  MatrixXd mwdcoscoefs,mwdsincoefs;
  MWD->select_coefficient_time(0.0, mwdcoscoefs, mwdsincoefs);

  double rtmp,phitmp,thetatmp,r2tmp;
  double tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  // the density fields are zombies right now: to be fixed with proper demand
  double dens0 = 0;
  double denstmp = 0;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys,pphys,dphys;

  // compute spherical coordinates in the frame of the MW expansion
  cartesian_to_cylindrical(xvir, yvir, r2tmp, phitmp);

  // same procedure for the disc
  MWD->determine_fields_at_point_cyl(mwdcoscoefs,mwdsincoefs,
                                     r2tmp,phitmp,zvir,
                                     tpotl0,tpotl,
                                     fr,fp,fztmp,mwdharmonicflag);

  cylindrical_forces_to_cartesian(rtmp, phitmp,
                                  fr, fp,
                                  fxtmp, fytmp);

  // translate to physical units
  virial_to_physical_density(denstmp,dphys);
  virial_to_physical_potential(tpotl,pphys);
  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // return MW force

  // initialise output array
  std::vector<double> output(5);

  output[0] = fxphys;
  output[1] = fyphys;
  output[2] = fzphys;
  output[3] = dphys;
  output[4] = pphys;

  return output;

}


std::vector<double>  MWLMC::all_forces(double t, double x, double y, double z,
                                       bool globalframe,
                                       int mwhharmonicflag, int mwdharmonicflag, int lmcharmonicflag,
                                       bool verbose)
{
  /*

   */

  // zero out forces
  double fx = 0;
  double fy = 0;
  double fz = 0;

  std::vector<double> output(3);


  // translate all times and positions into exp virial units
  double tvir,xvir,yvir,zvir;
  physical_to_virial_time(t,tvir);
  physical_to_virial_length(x,y,z, xvir,yvir,zvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  // get the initial coefficient values: the time here is in tvir units, so always start with 0
  MatrixXd tcoefsmw,tcoefslmc;
  MW->select_coefficient_time(tvir, tcoefsmw);
  LMC->select_coefficient_time(tvir, tcoefslmc);

  MatrixXd mwdcoscoefs,mwdsincoefs;
  MWD->select_coefficient_time(tvir, mwdcoscoefs, mwdsincoefs);

  // initialise the centre vectors
  vector<double> zerocoords(3),mw_centre(3),lmc_centre(3),mwd_centre(3);

  // get the present-day MWD coordinates: the zero of the system
  return_centre(reference_time, MWD->orient, zerocoords);

  //cout << "Coordinate zero:" << setw(14) << zerocoords[0] << setw(14) << zerocoords[1] << setw(14) << zerocoords[2] << endl;

  // get the centres of the expansions at the specified times in exp reference space
  return_centre(tvir,  MW->orient,  mw_centre);
  return_centre(tvir, LMC->orient, lmc_centre);
  return_centre(tvir, MWD->orient, mwd_centre);

  // shift the expansion centres to the pericentre coordinate system
  for (int j=0;j<=2;j++) {
    mw_centre[j]  -= zerocoords[j];
    lmc_centre[j] -= zerocoords[j];
    mwd_centre[j] -= zerocoords[j];
  }

  if (verbose) {
    cout << "MW virial centre (x,y,z)=(" << mw_centre[0] << ","<< mw_centre[1] << ","<< mw_centre[2] << ")" <<endl;
    cout << "MWD virial centre (x,y,z)=(" << mwd_centre[0] << ","<< mwd_centre[1] << ","<< mwd_centre[2] << ")" <<endl;
    cout << "LMC virial centre (x,y,z)=(" << lmc_centre[0] << ","<< lmc_centre[1] << ","<< lmc_centre[2] << ")" <<endl;
  }

  double rtmp,phitmp,thetatmp,r2tmp;
  double tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys;

  // compute spherical coordinates in the frame of the MW expansion
  cartesian_to_spherical(xvir-mwd_centre[0], yvir-mwd_centre[1], zvir-mwd_centre[2], rtmp, phitmp, thetatmp);

  //cout << setw(14) << rtmp << setw(14) << phitmp << setw(14) << thetatmp << endl;

  MW->determine_fields_at_point_sph(tcoefsmw,
                                    rtmp,thetatmp,phitmp,
                                    tpotl0,tpotl,
                                    fr,ft,fp,mwhharmonicflag);

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  r2tmp = sqrt((xvir-mwd_centre[0])*(xvir-mwd_centre[0]) + (yvir-mwd_centre[1])*(yvir-mwd_centre[1]));

  // same procedure for the disc
  MWD->determine_fields_at_point_cyl(mwdcoscoefs,mwdsincoefs,
                                     r2tmp,phitmp,zvir-mwd_centre[2],
                                     tpotl0,tpotl,
                                     fr,fp,fztmp,mwdharmonicflag);

  cylindrical_forces_to_cartesian(rtmp, phitmp,
                                  fr, fp,
                                  fxtmp, fytmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  // same procedure for LMC
  cartesian_to_spherical(xvir-lmc_centre[0], yvir-lmc_centre[1], zvir-lmc_centre[2], rtmp, phitmp, thetatmp);

  //cout << setw(14) << rtmp << setw(14) << phitmp << setw(14) << thetatmp << endl;

  LMC->determine_fields_at_point_sph(tcoefslmc,
                                     rtmp,thetatmp,phitmp,
                                     tpotl0,tpotl,
                                     fr,ft,fp,lmcharmonicflag);

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  // reset to physical units
  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add LMC force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  output[0] = fx;
  output[1] = fy;
  output[2] = fz;

  return output;

}

MatrixXd MWLMC::get_trajectories(double dt, bool virial)
{
  // get the trajectories with the specified time sampling (returns native time sampling if not specified)
  int nint = reference_time/dt;

  MatrixXd trajectories;
  trajectories.resize(nint,13);

  double tphys,xphys;


  for (int n=0;n<nint;n++)
  {
    // always start in native time, as we know the simulation bounds
    std::vector<double> trajectorytmp = get_expansion_centres_virial(n*dt);

    if (virial) {

      trajectories(n,0) = n*dt;
      for (int j=0;j<12;j++) trajectories(n,j+1) = trajectorytmp[j];

    } else {

      // translate back to physical units
      virial_to_physical_time(n*dt,tphys);
      trajectories(n,0) = tphys;

      for (int j=0;j<12;j++) {
        trajectories(n,j+1) = virial_to_physical_length(trajectorytmp[j]);
      }

    }
  }

  return trajectories;

}

MatrixXd MWLMC::get_lmc_trajectory(double dt)
{
  // get the trajectories with the specified time sampling (returns native time sampling if not specified)
  int nint = reference_time/dt;

  MatrixXd trajectory;
  trajectory.resize(nint,4);

  double tphys,xphys;

  for (int n=0;n<nint;n++)
  {
    // always start in native time, as we know the simulation bounds
    std::vector<double> trajectorytmp = get_lmc_centre_virial(n*dt,false);
    // translate back to physical units
    virial_to_physical_time(n*dt,tphys);
    trajectory(n,0) = tphys;

    for (int j=0;j<3;j++) {
      trajectory(n,j+1) = virial_to_physical_length(trajectorytmp[j]);
    }
  }

  return trajectory;

}


std::vector<double> MWLMC::get_expansion_centres_physical(double t, bool verbose)
{
  // check for a valid time
  if (t > 0.0) {
    std::cout << "Cannot select a time after the present day! Setting to present day..." << std::endl;
    t = 0.0;
  }

  // make one big return container
  std::vector<double> centres(12);

  // translate all times and positions into exp virial units
  double tvir;
  physical_to_virial_time(t,tvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  // initialise temporary centre vectors
  vector<double> zerocoords(3),mw_centre(3),lmc_centre(3),mwd_centre(3);

  // get the present-day MWD coordinates: the zero of the total
  return_centre(tvir, MWD->orient, zerocoords);

  if (verbose) std::cout << "Coordinate zero:" << setw(14) << zerocoords[0] << setw(14) << zerocoords[1] << setw(14) << zerocoords[2] << std::endl;

  // get the centres of the expansions at the specified times in exp reference space
  return_centre(tvir,  MW->orient,  mw_centre);
  return_centre(tvir, LMC->orient, lmc_centre);
  return_centre(tvir, MWD->orient, mwd_centre);

  // fill return vector
  for (int i=0;i<3;i++) centres[i  ] = zerocoords[i];
  for (int i=0;i<3;i++) centres[i+3] = mw_centre[i];
  for (int i=0;i<3;i++) centres[i+6] = lmc_centre[i];
  for (int i=0;i<3;i++) centres[i+9] = mwd_centre[i];

  return centres;

}


std::vector<double> MWLMC::get_expansion_centres_virial(double tvir, bool verbose)
{
  // check for a valid time
  if (tvir > reference_time) {
    std::cout << "Cannot select a time after the present day! Setting to present day..." << std::endl;
    tvir = reference_time;
  }

  // make one big return container
  std::vector<double> centres(12);

  // initialise temporary centre vectors
  vector<double> zerocoords(3),mw_centre(3),lmc_centre(3),mwd_centre(3);

  // get the present-day MWD coordinates: the zero of the total
  return_centre(tvir, MWD->orient, zerocoords);

  if (verbose) std::cout << "Coordinate zero:" << setw(14) << zerocoords[0] << setw(14) << zerocoords[1] << setw(14) << zerocoords[2] << std::endl;

  // get the centres of the expansions at the specified times in exp reference space
  return_centre(tvir,  MW->orient,  mw_centre);
  return_centre(tvir, LMC->orient, lmc_centre);
  return_centre(tvir, MWD->orient, mwd_centre);

  int indx;
  double dt;
  find_time_index(tvir, LMC->orient, indx, dt);

  // fill return vector
  for (int i=0;i<3;i++) centres[i  ] = zerocoords[i];
  for (int i=0;i<3;i++) centres[i+3] = mw_centre[i];
  for (int i=0;i<3;i++) centres[i+6] = lmc_centre[i];
  for (int i=0;i<3;i++) centres[i+9] = mwd_centre[i];

  return centres;

}


std::vector<double> MWLMC::get_lmc_centre_virial(double tvir, bool verbose)
{
  // check for a valid time
  if (tvir > reference_time) {
    std::cout << "Cannot select a time after the present day! Setting to present day..." << std::endl;
    tvir = reference_time;
  }

  // make one big return container
  std::vector<double> centres(3);

  // initialise temporary centre vectors
  vector<double> zerocoords(3),lmc_centre(3);

  // get the present-day MWD coordinates: the zero of the total
  return_centre(tvir, MWD->orient, zerocoords);

  // get the centres of the expansions at the specified times in exp reference space
  return_centre(tvir, LMC->orient, lmc_centre);

  // fill return vector
  for (int i=0;i<3;i++) centres[i  ] = lmc_centre[i] - zerocoords[i];

  return centres;

}


void MWLMC::all_forces_coefs(MatrixXd mwcoefs, MatrixXd lmccoefs,
                             MatrixXd mwdcoscoefs, MatrixXd mwdsincoefs,
                             double t, double x, double y, double z,
                             double& fx, double& fy, double& fz,
                             int mwhharmonicflag, int mwdharmonicflag, int lmcharmonicflag,
                             bool verbose)
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

  // initialise the centre vectors
  vector<double> zerocoords(3),mw_centre(3),lmc_centre(3),mwd_centre(3);
  return_centre(reference_time, MWD->orient, zerocoords); // the zero of the coordinate system
  return_centre(tvir,  MW->orient,  mw_centre);
  return_centre(tvir, LMC->orient, lmc_centre);
  return_centre(tvir, MWD->orient, mwd_centre);

  // shift the expansion centres to the pericentre coordinate system
  for (int j=0;j<=2;j++) {
    mw_centre[j]  -= zerocoords[j];
    lmc_centre[j] -= zerocoords[j];
    mwd_centre[j] -= zerocoords[j];
  }

  if (verbose) {
    cout << "MW virial centre (x,y,z)=(" << mw_centre[0] << ","<< mw_centre[1] << ","<< mw_centre[2] << ")" <<endl;
    cout << "MWD virial centre (x,y,z)=(" << mwd_centre[0] << ","<< mwd_centre[1] << ","<< mwd_centre[2] << ")" <<endl;
    cout << "LMC virial centre (x,y,z)=(" << lmc_centre[0] << ","<< lmc_centre[1] << ","<< lmc_centre[2] << ")" <<endl;
  }

  double rtmp,phitmp,thetatmp,r2tmp;
  double tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys;

  // compute spherical coordinates in the frame of the MW expansion
  cartesian_to_spherical(xvir-mwd_centre[0], yvir-mwd_centre[1], zvir-mwd_centre[2], rtmp, phitmp, thetatmp);

  //cout << setw(14) << rtmp << setw(14) << phitmp << setw(14) << thetatmp << endl;

  MW->determine_fields_at_point_sph(mwcoefs,
                                    rtmp,thetatmp,phitmp,
                                    tpotl0,tpotl,
                                    fr,ft,fp,mwhharmonicflag);

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  r2tmp = sqrt((xvir-mwd_centre[0])*(xvir-mwd_centre[0]) + (yvir-mwd_centre[1])*(yvir-mwd_centre[1]));

  // same procedure for the disc
  MWD->determine_fields_at_point_cyl(mwdcoscoefs,mwdsincoefs,
                                     r2tmp,phitmp,zvir-mwd_centre[2],
                                     tpotl0,tpotl,
                                     fr,fp,fztmp,mwdharmonicflag);

  cylindrical_forces_to_cartesian(rtmp, phitmp,
                                  fr, fp,
                                  fxtmp, fytmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  // same procedure for LMC
  cartesian_to_spherical(xvir-lmc_centre[0], yvir-lmc_centre[1], zvir-lmc_centre[2], rtmp, phitmp, thetatmp);

  //cout << setw(14) << rtmp << setw(14) << phitmp << setw(14) << thetatmp << endl;

  LMC->determine_fields_at_point_sph(lmccoefs,
                                     rtmp,thetatmp,phitmp,
                                     tpotl0,tpotl,
                                     fr,ft,fp,lmcharmonicflag);

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


void MWLMC::mw_forces_coefs(MatrixXd mwcoefs,
                            MatrixXd mwdcoscoefs, MatrixXd mwdsincoefs,
                            double t, double x, double y, double z,
                            double& fx, double& fy, double& fz,
                            int mwhharmonicflag, int mwdharmonicflag,
                            bool verbose)
{
  /*
    No translation: will only ever happen in the local frame of the MW.

   */

  // zero out forces
  fx = fy = fz = 0.;

  // translate all times and positions into exp virial units
  double xvir,yvir,zvir;
  physical_to_virial_length(x,y,z, xvir,yvir,zvir);

  double rtmp,phitmp,thetatmp,r2tmp;
  double tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys;

  // compute spherical coordinates
  cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);

  MW->determine_fields_at_point_sph(mwcoefs,
                                    rtmp,thetatmp,phitmp,
                                    tpotl0,tpotl,
                                    fr,ft,fp,mwhharmonicflag);

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  r2tmp = sqrt(xvir*xvir + yvir*yvir);

  // same procedure for the disc
  MWD->determine_fields_at_point_cyl(mwdcoscoefs,mwdsincoefs,
                                     r2tmp,phitmp,zvir,
                                     tpotl0,tpotl,
                                     fr,fp,fztmp,mwdharmonicflag);

  cylindrical_forces_to_cartesian(rtmp, phitmp,
                                  fr, fp,
                                  fxtmp, fytmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

}



MatrixXd MWLMC::orbit(vector<double> xinit,
                      vector<double> vinit,
                      double tbegin,
                      double tend,
                      double dt,
                      int mwhharmonicflag, int mwdharmonicflag, int lmcharmonicflag,
                      bool discframe)
{
  /*
  takes physical units

   */
   // start by converting everything to system units

  double fx,fy,fz,tvir;

  MatrixXd orbit;

  // number of integration steps to take
  int nint = (int)((tend - tbegin)/dt);

  //std::cout << "nint=" << nint << std::endl;

  if (nint<0) {
    std::cerr << "orbit: ending time must be before beginning time. Exiting." << std::endl;
    exit(-1);
  }

  // include the forces for now
  orbit.resize(10,nint);


  // set times
  double tvirbegin  = physical_to_virial_time_return(tbegin) + reference_time;
  double dtvir      = physical_to_virial_time_return(dt);
  double tphysbegin = tbegin;
  double dtphys     = dt;

  // initialise beginning values
  /*
  if (discframe) {
    vector<double> zerocoords(3),zerovel(3);
    // get the present-day MWD coordinates: the zero of the total
    return_centre(tvirbegin, MWD->orient, zerocoords);
    return_vel_centre(tvirbegin, MWD->orient, zerovel);
    for (j=0; j<3; j++) {
      orbit(j,n)   = orbit(j,n) - virial_to_physical_length(zerocoords[j]);
      orbit(j+3,n) = orbit(j,n) - virial_to_physical_velocity(zerovel[j]);
  }
  */


  orbit(0,0) = xinit[0];
  orbit(1,0) = xinit[1];
  orbit(2,0) = xinit[2];
  orbit(3,0) = vinit[0];
  orbit(4,0) = vinit[1];
  orbit(5,0) = vinit[2];
  //now step forward one, using leapfrog (drift-kick-drift) integrator?
  //    https://en.wikipedia.org/wiki/Leapfrog_integration
  //
  int step = 1;


  if (tvirbegin < 0) {
    std::cout << "fullintegrate.h: tvirbegin is less than 0! --> " << tvirbegin << std::endl;
  }

  // get the initial coefficient values: the time here is in tvir units, so always start with 0
  MatrixXd tcoefsmw,tcoefslmc, mwcoscoefs,mwsincoefs;
  MW->select_coefficient_time(tvirbegin, tcoefsmw);
  LMC->select_coefficient_time(tvirbegin, tcoefslmc);
  MWD->select_coefficient_time(tvirbegin, mwcoscoefs, mwsincoefs);

  // return forces for the initial step
  all_forces_coefs(tcoefsmw, tcoefslmc, mwcoscoefs, mwsincoefs,
                   tphysbegin, orbit(0,0),orbit(1,0),orbit(2,0),
                   fx, fy, fz,
                   mwhharmonicflag, mwdharmonicflag, lmcharmonicflag);

  orbit(6,0) = fx;
  orbit(7,0) = fy;
  orbit(8,0) = fz;

  int j;
  double tvirnow, tphysnow;

  for (step=1; step<nint; step++) {

    // advance timestep: this is in virial units by definition.
    tvirnow  = tvirbegin  + dtvir*step;
    tphysnow = tphysbegin + dtphys*step;
    orbit(9,step) = tphysnow; // record the time

    //std::cout << "step=" << step << " tphysnow=" << tphysnow << " tvirnow=" << tvirnow;

    MW->select_coefficient_time(tvirnow, tcoefsmw);
    LMC->select_coefficient_time(tvirnow, tcoefslmc);
    MWD->select_coefficient_time(tvirnow, mwcoscoefs, mwsincoefs);

    //std::cout << "...got coefficients";

    // advance positions
    for (j=0; j<3; j++) {
      orbit(j,step) = orbit(j,step-1)   + (orbit(j+3,step-1)*dt  )  + (0.5*orbit(j+6,step-1)  * (dt*dt));
    }

    // this call goes out in physical units
    all_forces_coefs(tcoefsmw, tcoefslmc, mwcoscoefs, mwsincoefs,
                     // need to shift the time by one unit to get the proper centre
                     tphysnow - dtphys, orbit(0,step),orbit(1,step),orbit(2,step),
                     fx, fy, fz,
                     mwhharmonicflag, mwdharmonicflag, lmcharmonicflag);

    //std::cout << "...got forces";


    orbit(6,step) = fx;
    orbit(7,step) = fy;
    orbit(8,step) = fz;

    // advance velocities
    for (j=3; j<6; j++) {
      orbit(j,step) = orbit(j,step-1) + (0.5*(orbit(j+3,step-1)+orbit(j+3,step))  * dt );
    }

    //std::cout << "...advanced velocities" << std::endl;

  }


  // convert to physical units again...
  if (discframe) {
    vector<double> zerocoords(3),zerovel(3);
    for (int n=0; n<nint; n++) {
      tvirnow  = tvirbegin  + dtvir*step;
      // get the present-day MWD coordinates: the zero of the total
      return_centre(tvirnow, MWD->orient, zerocoords);
      return_vel_centre(tvirnow, MWD->orient, zerovel);
      for (j=0; j<3; j++) {
        orbit(j,n)   = orbit(j,n) - virial_to_physical_length(zerocoords[j]);
        orbit(j+3,n) = orbit(j,n) - virial_to_physical_velocity(zerovel[j]);
        //orbit(j,n) = virial_to_physical_length(orbit(j,n));
        //orbit(j+3,n) = virial_to_physical_velocity(orbit(j+3,n));
        //orbit(j+6,n) = virial_to_physical_force(orbit(j+6,n));
      }
      //orbit(9,n) = virial_to_physical_time_return(orbit(9,n) - reference_time);
    }
  }

  return orbit;

}


MatrixXd MWLMC::mworbit(vector<double> xinit,
                        vector<double> vinit,
                        double tbegin,
                        double tend,
                        double dt)
{
  /*
  integration of just the Milky Way components
   */

   double fx,fy,fz,tvir;

   MatrixXd orbit;

   // number of integration steps to take
   int nint = (int)((tend - tbegin)/dt);


   if (nint<0) {
     std::cerr << "orbit: ending time must be before beginning time. Exiting." << std::endl;
     exit(-1);
   }

   // include the forces for now
   orbit.resize(10,nint);

   // set times
   double tvirbegin  = physical_to_virial_time_return(tbegin) + reference_time;
   double dtvir      = physical_to_virial_time_return(dt);
   double tphysbegin = tbegin;
   double dtphys     = dt;

   orbit(0,0) = xinit[0];
   orbit(1,0) = xinit[1];
   orbit(2,0) = xinit[2];
   orbit(3,0) = vinit[0];
   orbit(4,0) = vinit[1];
   orbit(5,0) = vinit[2];
   //now step forward one, using leapfrog (drift-kick-drift) integrator?
   //    https://en.wikipedia.org/wiki/Leapfrog_integration
   //
   int step = 1;


   if (tvirbegin < 0) {
     std::cout << "fullintegrate.h: tvirbegin is less than 0! --> " << tvirbegin << std::endl;
   }

   // get the initial coefficient values: the time here is in tvir units, so always start with 0
   MatrixXd tcoefsmw,tcoefslmc, mwcoscoefs,mwsincoefs;
   MW->select_coefficient_time(tvirbegin, tcoefsmw);
   MWD->select_coefficient_time(tvirbegin, mwcoscoefs, mwsincoefs);

   // return forces for the initial step
   mw_forces_coefs(tcoefsmw, mwcoscoefs, mwsincoefs,
                    tphysbegin, orbit(0,0),orbit(1,0),orbit(2,0),
                    fx, fy, fz,
                    127,127);

   orbit(6,0) = fx;
   orbit(7,0) = fy;
   orbit(8,0) = fz;

   int j;
   double tvirnow, tphysnow;

   for (step=1; step<nint; step++) {

     // advance timestep: this is in virial units by definition.
     tvirnow  = tvirbegin  + dtvir*step;
     tphysnow = tphysbegin + dtphys*step;
     orbit(9,step) = tphysnow; // record the time

     // advance positions
     for (j=0; j<3; j++) {
       orbit(j,step) = orbit(j,step-1)   + (orbit(j+3,step-1)*dt  )  + (0.5*orbit(j+6,step-1)  * (dt*dt));
     }

     // this call goes out in physical units
     mw_forces_coefs(tcoefsmw, mwcoscoefs, mwsincoefs,
                      // need to shift the time by one unit to get the proper centre
                      tphysbegin, orbit(0,step),orbit(1,step),orbit(2,step),
                      fx, fy, fz,
                      127,127);

     orbit(6,step) = fx;
     orbit(7,step) = fy;
     orbit(8,step) = fz;

     // advance velocities
     for (j=3; j<6; j++) {
       orbit(j,step) = orbit(j,step-1) + (0.5*(orbit(j+3,step-1)+orbit(j+3,step))  * dt );
     }

   }

   return orbit;

 }





MatrixXd MWLMC::rewind(vector<double> xinit,
                       vector<double> vinit,
                       double dt,
                       int mwhharmonicflag, int mwdharmonicflag, int lmcharmonicflag,
                       double rewind_time,
                       bool discframe)
{
 /*
 takes physical units

  */
  // start by converting everything to system units

 double fx,fy,fz,tvir;

 MatrixXd orbit;

 // number of integration steps to take
 int nint = (int)(rewind_time/dt);

 if (nint<0) {
   std::cerr << "rewind: rewind_time must be positive. Exiting." << std::endl;
   exit(-1);
 }

 // include the forces for now
 orbit.resize(10,nint);

 // set times
 double tvirbegin  = reference_time;
 double dtvir      = physical_to_virial_time_return(dt);
 double tphysbegin = 0.;
 double dtphys     = dt;

 orbit(0,0) =  xinit[0];
 orbit(1,0) =  xinit[1];
 orbit(2,0) =  xinit[2];
 orbit(3,0) = -vinit[0];
 orbit(4,0) = -vinit[1];
 orbit(5,0) = -vinit[2];
 //now step forward one, using leapfrog (drift-kick-drift) integrator?
 //    https://en.wikipedia.org/wiki/Leapfrog_integration
 //
 int step = 1;

 // get the initial coefficient values: the time here is in tvir units, so always start with 0
 MatrixXd tcoefsmw,tcoefslmc, mwcoscoefs,mwsincoefs;
 MW->select_coefficient_time(tvirbegin, tcoefsmw);
 LMC->select_coefficient_time(tvirbegin, tcoefslmc);
 MWD->select_coefficient_time(tvirbegin, mwcoscoefs, mwsincoefs);

 // return forces for the initial step
 all_forces_coefs(tcoefsmw, tcoefslmc, mwcoscoefs, mwsincoefs,
                  tphysbegin, orbit(0,0),orbit(1,0),orbit(2,0),
                  fx, fy, fz,
                  mwhharmonicflag, mwdharmonicflag, lmcharmonicflag);

 orbit(6,0) = fx;
 orbit(7,0) = fy;
 orbit(8,0) = fz;

 int j;
 double tvirnow, tphysnow;

 for (step=1; step<nint; step++) {

   // decrement timestep: this is in virial units by definition.
   tvirnow  = tvirbegin  - dtvir*step;
   tphysnow = tphysbegin - dtphys*step;
   orbit(9,step) = tphysnow; // record the time

   //std::cout << "step=" << step << " tphysnow=" << tphysnow << " tvirnow=" << tvirnow << std::endl;

   MW->select_coefficient_time(tvirnow, tcoefsmw);
   LMC->select_coefficient_time(tvirnow, tcoefslmc);
   MWD->select_coefficient_time(tvirnow, mwcoscoefs, mwsincoefs);

   // 'advance' positions
   for (j=0; j<3; j++) {
     orbit(j,step) = orbit(j,step-1)   + (orbit(j+3,step-1)*dt  )  + (0.5*orbit(j+6,step-1)  * (dt*dt));
   }

   // this call goes out in physical units
   all_forces_coefs(tcoefsmw, tcoefslmc, mwcoscoefs, mwsincoefs,
                    // need to shift the time by one unit to get the proper centre
                    tphysnow - dtphys, orbit(0,step),orbit(1,step),orbit(2,step),
                    fx, fy, fz,
                    mwhharmonicflag, mwdharmonicflag, lmcharmonicflag);


   orbit(6,step) = fx;
   orbit(7,step) = fy;
   orbit(8,step) = fz;

   // advance velocities
   for (j=3; j<6; j++) {
     orbit(j,step) = orbit(j,step-1) + (0.5*(orbit(j+3,step-1)+orbit(j+3,step))  * dt );
   }

 }


 // convert to physical units again...
 vector<double> initcoords(3),initvel(3);
 tvirnow  = tvirbegin;
 return_centre(tvirnow, MWD->orient, initcoords);
 //std::cout << setw(14) << tvirnow  << setw(14) << zerocoords[0] << setw(14) << zerocoords[1] << setw(14) << zerocoords[2] << std::endl;
 return_vel_centre(tvirnow, MWD->orient, initvel);

 if (discframe) {
   vector<double> zerocoords(3),zerovel(3);
   for (int n=0; n<nint; n++) {
     tvirnow  = tvirbegin - dtvir*n;
     // get the present-day MWD coordinates: the zero of the total
     return_centre(tvirnow, MWD->orient, zerocoords);
     //std::cout << setw(14) << tvirnow  << setw(14) << zerocoords[0] << setw(14) << zerocoords[1] << setw(14) << zerocoords[2] << std::endl;
     return_vel_centre(tvirnow, MWD->orient, zerovel);
     for (j=0; j<3; j++) {
       orbit(j,n)   = orbit(j,n) - virial_to_physical_length(zerocoords[j]) + virial_to_physical_length(initcoords[j]);
       orbit(j+3,n) = orbit(j,n) - virial_to_physical_velocity(zerovel[j]) + virial_to_physical_length(initvel[j]);
       //orbit(j,n) = virial_to_physical_length(orbit(j,n));
       //orbit(j+3,n) = virial_to_physical_velocity(orbit(j+3,n));
       //orbit(j+6,n) = virial_to_physical_force(orbit(j+6,n));
     }
     //orbit(9,n) = virial_to_physical_time_return(orbit(9,n) - reference_time);
   }
 }

 return orbit;

}


#endif
