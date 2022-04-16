/*
example code:
integrate an orbit in the MW-LMC potential and compare to exp output

compile string:
clang++ -I/opt/local/include -L/opt/local/lib -Iinclude/ fullintegrate.cc -o fullintegrate
clang++ -I/opt/local/include -L/opt/local/lib -I../include/ fullintegrate.cc -o obj/fullintegrate

clang++ --std=c++17 -lyaml-cpp -I/opt/local/include -L/opt/local/lib -I../include/ fullintegrate.cc -o obj/fullintegrate
otool -L obj/fullintegrate
install_name_tool -change @rpath/libyaml-cpp.0.7.dylib /opt/local/lib/libyaml-cpp.dylib obj/fullintegrate

clang++ --std=c++17 -lyaml-cpp -I/opt/local/include  -I../include/ fullintegrate.cc -o obj/fullintegrate


-O3 doesn't add much.
-g -O0 really slows down.

MSP 28 Apr 2020 initial commit
MSP 13 Oct 2020 new model validation
MSP 22 Dec 2021 add tests for new yaml formats

*/

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



void print_orbit(MatrixXd orbit,
		 string orbitfile)
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



void return_forces_mw_and_lmc(SphExpansion* MW, SphExpansion* LMC,
			      MatrixXd mwcoefs, MatrixXd lmccoefs,
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

  MW->determine_fields_at_point_sph(mwcoefs,
				    rtmp,thetatmp,phitmp,
				    tpotl0,tpotl,
				    fr,ft,fp,false,false,false,2);

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

  LMC->determine_fields_at_point_sph(lmccoefs,
				     rtmp,thetatmp,phitmp,
				     tpotl0,tpotl,
				     fr,ft,fp,true,true,true);

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
			    SphExpansion* LMC,
			    vector<double> xinit,
			    vector<double> vinit,
			    int nint,
			    double dt,
			    MatrixXd& orbit)
{
  /*

   */
  double fx,fy,fz;


  double tvir;

  // include the forces for now
  orbit.resize(10,nint);

	cout << "Orbit dims " << orbit.rows() << " x " << orbit.cols() << endl;

  // initialise beginning values
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

  // get the initial coefficient values: the time here is in tvir units, so always start with 0
  MatrixXd tcoefsmw,tcoefslmc;
  MW->select_coefficient_time(0.,tcoefsmw);
  LMC->select_coefficient_time(0., tcoefslmc);

  // not applying time offsets here; think about whether this is a problem
  double tphys;
  virial_to_physical_time(0.,tphys);
  // return forces for the initial step
  return_forces_mw_and_lmc(MW, LMC,
			     tcoefsmw, tcoefslmc,
			     tphys, orbit(0,0),orbit(1,0),orbit(2,0),
			     fx, fy, fz, false);

  orbit(6,0) = fx;
  orbit(7,0) = fy;
  orbit(8,0) = fz;

  int j;

  for (step=1; step<nint; step++) {

    // advance timestep: this is in physical units by definition.
    orbit(9,step) = dt*step;

    // find the current virial time
    physical_to_virial_time(dt*(step),tvir);

    // check that this is incrementing correctly? it is.
    //cout << "TVIR=" << tvir << endl;

    // get coefficients at the current virial time
    MW->select_coefficient_time(tvir, tcoefsmw);
    LMC->select_coefficient_time(tvir,  tcoefslmc);

    // advance positions
    for (j=0; j<3; j++) {
      orbit(j,step) = orbit(j,step-1)   + (orbit(j+3,step-1)*dt  )  + (0.5*orbit(j+6,step-1)  * (dt*dt));
    }

    // calculate new forces: time goes in as physical time (e.g. kpc/km/s)
    return_forces_mw_and_lmc(MW, LMC,
			     tcoefsmw, tcoefslmc,
			     dt*(step-1), orbit(0,step),orbit(1,step),orbit(2,step),
			     fx, fy, fz, false);

    orbit(6,step) = fx;
    orbit(7,step) = fy;
    orbit(8,step) = fz;

    // advance velocities
    for (j=3; j<6; j++) {
      orbit(j,step) = orbit(j,step-1) + (0.5*(orbit(j+3,step-1)+orbit(j+3,step))  * dt );
    }

  }

}



void return_forces_mw_and_lmc_with_disc(SphExpansion* MW, SphExpansion* LMC, CylExpansion* MWD,
			                MatrixXd mwcoefs, MatrixXd lmccoefs, MatrixXd mwdcoscoefs, MatrixXd mwdsincoefs,
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

  MW->determine_fields_at_point_sph( mwcoefs,
				rtmp,thetatmp,phitmp,
				tpotl0,tpotl,
				     fr,ft,fp,false,false,false);

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
				     fr,fp,fztmp,false,false,false);

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
				     fr,ft,fp,false,false,false);

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



void three_component_leapfrog(SphExpansion* MW,
			      SphExpansion* LMC,
			      CylExpansion* MWD,
			      vector<double> xinit,
			      vector<double> vinit,
			      int nint,
			      double dt,
			      MatrixXd& orbit)
{
  /*

   */
  double fx,fy,fz;


  double tvir;

  // include the forces for now
  orbit.resize(10,nint);

	cout << "Orbit dims " << orbit.rows() << " x " << orbit.cols() << endl;


  // initialise beginning values
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

  // get the initial coefficient values: the time here is in tvir units, so always start with 0
  MatrixXd tcoefsmw,tcoefslmc;
  MW->select_coefficient_time(0., tcoefsmw,15);
  LMC->select_coefficient_time(0., tcoefslmc);

  MatrixXd mwcoscoefs,mwsincoefs;
  MWD->select_coefficient_time(0.0, mwcoscoefs, mwsincoefs);

  // not applying time offsets here; think about whether this is a problem
  double tphys;
  virial_to_physical_time(0.,tphys);
  // return forces for the initial step
  return_forces_mw_and_lmc_with_disc(MW, LMC, MWD,
				     tcoefsmw, tcoefslmc, mwcoscoefs, mwsincoefs,
				     tphys, orbit(0,0),orbit(1,0),orbit(2,0),
				     fx, fy, fz, false);

  orbit(6,0) = fx;
  orbit(7,0) = fy;
  orbit(8,0) = fz;

  int j;

  for (step=1; step<nint; step++) {

    // advance timestep: this is in physical units by definition.
    orbit(9,step) = dt*step;

    // find the current virial time
    physical_to_virial_time(dt*(step),tvir);

    // check that this is incrementing correctly? it is.
    //cout << "TVIR=" << tvir << endl;

    // get coefficients at the current virial time
    MW->select_coefficient_time(tvir, tcoefsmw);
    LMC->select_coefficient_time(tvir, tcoefslmc);
    MWD->select_coefficient_time(tvir, mwcoscoefs, mwsincoefs);

    // advance positions
    for (j=0; j<3; j++) {
      orbit(j,step) = orbit(j,step-1)   + (orbit(j+3,step-1)*dt  )  + (0.5*orbit(j+6,step-1)  * (dt*dt));
    }

    // calculate new forces: time goes in as physical time (e.g. kpc/km/s)
    return_forces_mw_and_lmc_with_disc(MW, LMC, MWD,
				       tcoefsmw, tcoefslmc, mwcoscoefs, mwsincoefs,
				       dt*(step-1), orbit(0,step),orbit(1,step),orbit(2,step),
				       fx, fy, fz, false);

		orbit(6,step) = fx;
		orbit(7,step) = fy;
		orbit(8,step) = fz;

		// advance velocities
		for (j=3; j<6; j++) {
			orbit(j,step) = orbit(j,step-1) + (0.5*(orbit(j+3,step-1)+orbit(j+3,step))  * dt );
		}

  }

}






int main () {

  // MW
  cout << "Initialising MW ... " << endl;
  SphExpansion* MW = new SphExpansion(sph_cache_name_mw, model_file_mw, coef_file_mw, orient_file_mw);

  // LMC
  cout << "Initialising LMC ... " << endl;
  SphExpansion* LMC = new SphExpansion(sph_cache_name_lmc, model_file_lmc, coef_file_lmc, orient_file_lmc);

  // MW
  cout << "Initialising MW disc ... " << endl;
  CylExpansion* MWD = new CylExpansion(cyl_cache_name_mw, cyl_coef_name_mw, cyl_orient_name_mw);

  vector<double> zerocoords(3);
  return_centre(reference_time, MWD->orient, zerocoords);

  cout << "Coordinate zero:" << setw(14) << zerocoords[0] << setw(14) << zerocoords[1] << setw(14) << zerocoords[2] << endl;


  // get the time in exp system units by 1) converting the virial units, 2) applying time offset for present-day

  // initial position offset from centre
  vector<double> xphys(3),xinit(3);
  vector<double> vxphys(3),vxinit(3);

  // LMC ORBIT 1: reproduced
  //-0.0276161 1.77453 -0.0771715
  xphys[0] = -0.0276161;
  xphys[1] = 1.77453;
  xphys[2] = -0.0771715;
  //-0.0451152 -0.351639 -0.199325
  vxphys[0] = -0.0451152;
  vxphys[1] = -0.351639;
  vxphys[2] = -0.199325;

  // LMC ORBIT 2: Closer to the centre
  xphys[0] = 0.226053;
  xphys[1] = 2.09839;
  xphys[2] = -0.207248;
  vxphys[0] = 0.0510511;
  vxphys[1] = -0.458447;
  vxphys[2] = -0.24932;

  // MW ORBIT 1: Big loops
  xphys[0] = 0.123143;
  xphys[1] = 0.558696;
  xphys[2] = 0.0459567;
  vxphys[0] = 0.457418;
  vxphys[1] = -0.464346;
  vxphys[2] = -0.124515;

  // MW ORBIT z: Disc-like orbit
  xphys[0] = 0.03;
  xphys[1] = 0.0;
  xphys[2] = 0.0;
  vxphys[0] = 0.0;
  vxphys[1] = 1.4;
  vxphys[2] = 0.;

  virial_to_physical_length(xphys[0],xphys[1],xphys[2],xinit[0],xinit[1],xinit[2]);
  virial_to_physical_velocity(vxphys[0],vxphys[1],vxphys[2],vxinit[0],vxinit[1],vxinit[2]);

  cout << "Input pos/vel: " << xinit[0] << " " << xinit[1] << " " << xinit[2] << " " <<
    vxinit[0] << " " << vxinit[1] << " " << vxinit[2] << " " << endl;

  double nint=100;

  // call this time in kpc/km/s units
  double dt;
  // force sampling at the native rate of the exp simulation as an interpolation test
  virial_to_physical_time(0.0005,dt);
  MatrixXd orbit;

  //two_component_leapfrog(MW, LMC, xinit, vxinit, nint, dt, orbit);

  string orbitfile="tests/comparisonorbit6.txt";
  //print_orbit(orbit,orbitfile);

  // try with the disk as well!
  three_component_leapfrog(MW, LMC, MWD, xinit, vxinit, nint, dt, orbit);

  orbitfile="tests/comparisonorbit7.txt";
  //print_orbit(orbit,orbitfile);

}
