/*
tests for the initial MW and LMC
- density under deforming conditions
- made to track down a bug.

compile string:
clang++ --std=c++17 -I/opt/local/include -L/opt/local/lib -I../include/ density_map.cc -o obj/density_map

MSP 13 Oct 2020 first written.
MSP  4 Nov 2021 rewritten for new expansion.h
MSP 16 Apr 2022 point to stable example files

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

// expansion includes
#include "sphexpansion.h"


void make_surface(SphExpansion* S1,
		  MatrixXd coefs1,
		  SphExpansion* S2,
		  MatrixXd coefs2,
		  double reftime,
		  double xmin,
		  double xmax,
		  int nsamples,
		  string outfile,
		  int harmonicflag=SPHHARMONICDEFAULT)
{
  /*
  xmin and xmax come in PHYSICAL units
  and the density is returned in physical units
   */

  // hardwire the number of angular samples
  int asamples = nsamples;

  vector<double> mwcentre(3),lmccentre(3);
  return_centre(reftime,  S1->orient,  mwcentre);
  double mwxvp,mwyvp,mwzvp;
  virial_to_physical_length(mwcentre[0],mwcentre[1],mwcentre[2],mwxvp,mwyvp,mwzvp);

  return_centre(reftime,  S2->orient, lmccentre);
  double lmcxvp,lmcyvp,lmczvp;
  virial_to_physical_length(lmccentre[0],lmccentre[1],lmccentre[2],lmcxvp,lmcyvp,lmczvp);

  cout << setw(14) <<  mwxvp << setw(14) <<  mwyvp << setw(14) <<  mwzvp << endl;
  cout << setw(14) << lmcxvp << setw(14) << lmcyvp << setw(14) << lmczvp << endl;


  ofstream densityfile;
  densityfile.open(outfile);

  densityfile << "#  r [kpc] ; phi [deg] ; theta [deg] ;mwdensity [Msun/pc^3] ; lmcdensity [Msun/pc^3];" << endl;

  double dx = (xmax-xmin)/nsamples;
  double da = (2*3.14159)/asamples;
  double xin,yin,zin;
  double rin,tin,pin;
  double mwd,lmcd,mwphysdens,lmcphysdens;

  // empirically check the transformations
  spherical_to_cartesian(1.,0.,0.,
			         xin,yin,zin);

  cout << setw(14) << xin << setw(14) << yin  << setw(14) << zin << endl;

  spherical_to_cartesian(1.,0.,3.1415926,
			         xin,yin,zin);

  cout << setw(14) << xin << setw(14) << yin  << setw(14) << zin << endl;


  // demo to make a rotation curve
  for (int xx=0; xx<nsamples; xx++) {

    // the radius of the shell to probe
    rin = xx*dx + xmin; // in kpc

    for (int yy=0; yy<asamples; yy++) {

      pin = yy*da;

      for (int zz=0; zz<asamples; zz++) {

	  			tin = zz*0.5*da;

          spherical_to_cartesian(rin,pin,tin,
			         xin,yin,zin);


          S1->return_density(coefs1,
			     xin, yin, zin,
		             mwd, harmonicflag);

          S2->return_density(coefs2,
			     xin-lmcxvp, yin-lmcyvp, zin-lmczvp,
		             lmcd, harmonicflag);

          virial_to_physical_density( mwd,  mwphysdens);
          virial_to_physical_density(lmcd, lmcphysdens);

          densityfile  << setw(14) << rin        << setw(14) << pin << setw(14) << tin <<
                          setw(14) << mwphysdens << setw(14) << lmcphysdens << endl;

      } // end of z/theta loop
    } // end of y/phi loop
  } // end of x/r loop

  densityfile.close();

}


void make_density(SphExpansion* S1,
		  MatrixXd coefs1,
		  SphExpansion* S2,
		  MatrixXd coefs2,
		  double reftime,
		  double xmin,
		  double xmax,
		  int nsamples,
		  string outfile,
		  int harmonicflag=SPHHARMONICDEFAULT)
{
  /*
  xmin and xmax come in PHYSICAL units
  and the rotation curve is returned in physical units
   */

  vector<double> mwcentre(3),lmccentre(3);
  return_centre(reftime,  S1->orient,  mwcentre);
  double mwxvp,mwyvp,mwzvp;
  virial_to_physical_length(mwcentre[0],mwcentre[1],mwcentre[2],mwxvp,mwyvp,mwzvp);

  return_centre(reftime,  S2->orient, lmccentre);
  double lmcxvp,lmcyvp,lmczvp;
  virial_to_physical_length(lmccentre[0],lmccentre[1],lmccentre[2],lmcxvp,lmcyvp,lmczvp);

  cout << setw(14) <<  mwxvp << setw(14) <<  mwyvp << setw(14) <<  mwzvp << endl;
  cout << setw(14) << lmcxvp << setw(14) << lmcyvp << setw(14) << lmczvp << endl;


  ofstream densityfile;
  densityfile.open(outfile);

  densityfile << "#  y [kpc] ; z [km/s/s] ; mwdensity [Msun/pc^3] ; lmcdensity [Msun/pc^3];" << endl;

  double dx = (xmax-xmin)/nsamples;
  double zin,yin;
  double mwd,lmcd,mwphysdens,lmcphysdens;

  // demo to make a rotation curve
  for (int xx=0; xx<nsamples; xx++) {

    // the location in inertial space of the points to check (xin=0, just checking y-z plane right now)
    zin = xx*dx + xmin; // in kpc

    for (int yy=0; yy<nsamples; yy++) {

      yin = yy*dx + xmin;


      S1->return_density(coefs1,
			 //0.-mwxvp, yin-mwyvp, zin-mwzvp,
			 0., yin, zin,
		         mwd, harmonicflag);

      S2->return_density(coefs2,
			 0.-lmcxvp, yin-lmcyvp, zin-lmczvp,
			 //0.-mwxvp, yin-mwyvp, zin-mwzvp,
		         lmcd, harmonicflag);

      virial_to_physical_density( mwd,  mwphysdens);
      virial_to_physical_density(lmcd, lmcphysdens);

      densityfile  << setw(14) << yin        << setw(14) << zin <<
                      setw(14) << mwphysdens << setw(14) << lmcphysdens << endl;

    }
  }

  densityfile.close();

}


void make_density_3d(SphExpansion* S1,
		     MatrixXd coefs1,
		     SphExpansion* S2,
		     MatrixXd coefs2,
		     double reftime,
		     double xmin,
		     double xmax,
		     int nsamples,
		     string outfile,
		     int harmonicflag=SPHHARMONICDEFAULT)
{
  /*
  xmin and xmax come in PHYSICAL units
  and the rotation curve is returned in physical units
   */

  vector<double> mwcentre(3),lmccentre(3);
  return_centre(reftime,  S1->orient,  mwcentre);
  double mwxvp,mwyvp,mwzvp;
  virial_to_physical_length(mwcentre[0],mwcentre[1],mwcentre[2],mwxvp,mwyvp,mwzvp);

  return_centre(reftime,  S2->orient, lmccentre);
  double lmcxvp,lmcyvp,lmczvp;
  virial_to_physical_length(lmccentre[0],lmccentre[1],lmccentre[2],lmcxvp,lmcyvp,lmczvp);

  cout << "Inertial space locations:" << endl;
  cout << setw(14) <<  mwxvp << setw(14) <<  mwyvp << setw(14) <<  mwzvp << endl;
  cout << setw(14) << lmcxvp << setw(14) << lmcyvp << setw(14) << lmczvp << endl;


  ofstream densityfile;
  densityfile.open(outfile);

  densityfile << "#  y [kpc] ; z [km/s/s] ; mwdensity [Msun/pc^3] ; lmcdensity [Msun/pc^3];" << endl;

  double dx = (xmax-xmin)/nsamples;
  double xin,yin,zin;
  double mwd,lmcd,mwphysdens,lmcphysdens;

  // demo to make a rotation curve
  for (int xx=0; xx<nsamples; xx++) {

    xin = xx*dx + xmin; // in kpc

    for (int yy=0; yy<nsamples; yy++) {

      yin = yy*dx + xmin;

      for (int zz=0; zz<nsamples; zz++) {

				zin = zz*dx + xmin;


				S1->return_density(coefs1,
						   xin, yin, zin,
						   mwd, harmonicflag);

				S2->return_density(coefs2,
						   xin-lmcxvp, yin-lmcyvp, zin-lmczvp,
						   lmcd, harmonicflag);

				virial_to_physical_density( mwd,  mwphysdens);
				virial_to_physical_density(lmcd, lmcphysdens);

				densityfile  << setw(14) << xin
					           << setw(14) << yin
					           << setw(14) << zin
					           << setw(14) << mwphysdens
					           << setw(14) << lmcphysdens
					           << endl;
      } // zloop

    } // yloop

  }// xloop

  densityfile.close();

} // function loop



int main () {


  // MW
  cout << "Initialising MW ... " << endl;

  SphExpansion* MW;
  MW = new SphExpansion(sph_cache_name_mw, model_file_mw, coef_file_mw, orient_file_mw);

  // LMC
  cout << "Initialising LMC ... " << endl;

  SphExpansion* LMC;
  LMC = new SphExpansion(sph_cache_name_lmc, model_file_lmc, coef_file_lmc, orient_file_lmc);

	// the reference time for the simulation (present day) comes from modelfiles.h
  MatrixXd mwcoefs,lmccoefs;
  MW->select_coefficient_time(reference_time, mwcoefs);
  LMC->select_coefficient_time(reference_time, lmccoefs);

  string dfile="/Users/mpetersen/Desktop/MW_LMC_density.txt";
  //make_density(MW,mwcoefs,LMC,lmccoefs,reftime,-50,50.,100,dfile);

  string dfile3d="/Users/mpetersen/Desktop/MW_LMC_density_3d_3.txt";
  //make_density_3d(MW,mwcoefs,LMC,lmccoefs,reftime,-120,120.,35,dfile3d);

  string dfilesurf="/Users/mpetersen/Desktop/MW_LMC_density_3d_4.txt";
  make_surface(MW,mwcoefs,LMC,lmccoefs,reference_time,0.,100.,35,dfilesurf);

}
