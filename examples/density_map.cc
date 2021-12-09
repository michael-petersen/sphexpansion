/*
tests for the initial MW and LMC
- density under deforming conditions
- made to track down a bug.

compile string: 
clang++ -I/opt/local/include -L/opt/local/lib -I../include/ density_map.cc -o obj/density_map

MSP 13 Oct 2020 first written.
MSP  4 Nov 2021 rewritten for new expansion.h

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdio.h>

// boost includes
#include "boost/multi_array.hpp"

// expansion includes
#include "expansion.h"

// integration include
#include "leapfrog.h"

void make_surface(SphExpansion* S1,
		  array_type2 coefs1,
		  SphExpansion* S2,
		  array_type2 coefs2,
		  double reftime,
		  double xmin,
		  double xmax,
		  int nsamples,
		  string outfile,
		  bool monopole=false,
		  bool dipole=false,
		  bool quadrupole=false,
		  int ltrunc=10)
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
		             mwd, monopole, dipole, quadrupole, ltrunc);
      
          S2->return_density(coefs2,
			     xin-lmcxvp, yin-lmcyvp, zin-lmczvp,
		             lmcd, monopole, dipole, quadrupole, ltrunc);
    
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
		  array_type2 coefs1,
		  SphExpansion* S2,
		  array_type2 coefs2,
		  double reftime,
		  double xmin,
		  double xmax,
		  int nsamples,
		  string outfile,
		  bool monopole=false,
		  bool dipole=false,
		  bool quadrupole=false,
		  int ltrunc=10)
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
		         mwd, monopole, dipole, quadrupole, ltrunc);
      
      S2->return_density(coefs2,
			 0.-lmcxvp, yin-lmcyvp, zin-lmczvp,
			 //0.-mwxvp, yin-mwyvp, zin-mwzvp,
		         lmcd, monopole, dipole, quadrupole, ltrunc);
    
      virial_to_physical_density( mwd,  mwphysdens);
      virial_to_physical_density(lmcd, lmcphysdens);
    
      densityfile  << setw(14) << yin        << setw(14) << zin <<
                      setw(14) << mwphysdens << setw(14) << lmcphysdens << endl;

    }
  }
  
  densityfile.close();
  
}
			 

void make_density_3d(SphExpansion* S1,
		     array_type2 coefs1,
		     SphExpansion* S2,
		     array_type2 coefs2,
		     double reftime,
		     double xmin,
		     double xmax,
		     int nsamples,
		     string outfile,
		     bool monopole=false,
		     bool dipole=false,
		     bool quadrupole=false,
		     int ltrunc=10)
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
			   //0.-mwxvp, yin-mwyvp, zin-mwzvp,
			   xin, yin, zin,
			   mwd, monopole, dipole, quadrupole, ltrunc);
      
	S2->return_density(coefs2,
			   xin-lmcxvp, yin-lmcyvp, zin-lmczvp,
			   //0.-mwxvp, yin-mwyvp, zin-mwzvp,
			   lmcd, monopole, dipole, quadrupole, ltrunc);
    
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
  /*
  string sph_cache_name_mw  = "/Volumes/External1/Disk080/SLGridSph.cache.mw.system1_3";
  string model_file_mw      = "/Volumes/External1/Disk080/SLGridSph.NFW77";
  //string coef_file_mw       = "/Volumes/External1/Disk080/simpleoutcoef.nofac.mw.system1_3";
  string coef_file_mw       = "/Volumes/External1/Disk080/simpleoutcoef.nfac.mw.system1_3";
  string orient_file_mw     = "/Volumes/External1/Disk080/mw.orient.system1_3.smth";
  */
  string sph_cache_name_mw  = "/Volumes/External1/Disk080/SLGridSph.cache.mw.system2_3";
  string model_file_mw      = "/Volumes/External1/Disk080/SLGridSph.Hern77";
  string coef_file_mw       = "/Volumes/External1/Disk080/simpleoutcoef.nfac.mw.system2_3";
  string orient_file_mw     = "/Volumes/External1/Disk080/mw.orient.system2_3.smth";

  
  SphExpansion* MW;
  MW = new SphExpansion(sph_cache_name_mw, model_file_mw, coef_file_mw, orient_file_mw);

  // LMC
  cout << "Initialising LMC ... " << endl;
  /*
  string sph_cache_name_lmc  = "/Volumes/External1/Disk080/SLGridSph.cache.lmc.system1_3";
  string model_file_lmc      = "/Volumes/External1/Disk080/SLGridSph.NFW77L";
  //string coef_file_lmc       = "/Volumes/External1/Disk080/simpleoutcoef.nofac.lmc.system1_3";
  string coef_file_lmc       = "/Volumes/External1/Disk080/simpleoutcoef.nfac.lmc.system1_3";
  string orient_file_lmc     =
  "/Volumes/External1/Disk080/lmc.orient.system1_3.smth";
  */

  string sph_cache_name_lmc  = "/Volumes/External1/Disk080/SLGridSph.cache.lmc.system2_3";
  string model_file_lmc      = "/Volumes/External1/Disk080/SLGridSph.NFW77L";
  //string coef_file_lmc       = "/Volumes/External1/Disk080/simpleoutcoef.nofac.lmc.system1_3";
  string coef_file_lmc       = "/Volumes/External1/Disk080/simpleoutcoef.nfac.lmc.system2_3";
  string orient_file_lmc     = "/Volumes/External1/Disk080/lmc.orient.system2_3.smth";
  
  SphExpansion* LMC;
  LMC = new SphExpansion(sph_cache_name_lmc, model_file_lmc, coef_file_lmc, orient_file_lmc);
  

  /* 
  // MW
  cout << "Initialising MW ... " << endl;
  string sph_cache_name_mw  = "/Volumes/External1/Disk076/SLGridSph.mw.run9mld";
  string model_file_mw      = "/Volumes/External1/Disk076/ErkalMW.model";
  string coef_file_mw       = "/Volumes/External1/Disk076/simpleoutcoef.nofac.mw.run9mld";
  string orient_file_mw     = "/Volumes/External1/Disk076/mw.orient.run9mld.smth";

  SphExpansion* MW;
  MW = new SphExpansion(sph_cache_name_mw, model_file_mw, coef_file_mw, orient_file_mw);

  // LMC
  cout << "Initialising LMC ... " << endl;
  string sph_cache_name_lmc = "/Volumes/External1/Disk076/SLGridSph.lmc.run9mld";
  string model_file_lmc     = "/Volumes/External1/Disk076/ErkalLMC.model";
  string coef_file_lmc      = "/Volumes/External1/Disk076/simpleoutcoef.nofac.lmc.run9mld";
  string orient_file_lmc    = "/Volumes/External1/Disk076/lmc.orient.run9mld.smth";

  SphExpansion* LMC;
  LMC = new SphExpansion(sph_cache_name_lmc, model_file_lmc, coef_file_lmc, orient_file_lmc);
  */

  double reftime = 1.19;//1.393;
  reftime = 1.393;

  array_type2 mwcoefs,lmccoefs;
  MW->select_coefficient_time(reftime, mwcoefs);
  LMC->select_coefficient_time(reftime, lmccoefs);
  
  string dfile="/Users/mpetersen/Desktop/MW_LMC_density.txt";
  //make_density(MW,mwcoefs,LMC,lmccoefs,reftime,-50,50.,100,dfile);

  string dfile3d="/Users/mpetersen/Desktop/MW_LMC_density_3d_3.txt";
  //make_density_3d(MW,mwcoefs,LMC,lmccoefs,reftime,-120,120.,35,dfile3d);

  string dfilesurf="/Users/mpetersen/Desktop/MW_LMC_density_3d_4.txt";
  make_surface(MW,mwcoefs,LMC,lmccoefs,reftime,0.,100.,35,dfilesurf);
  
}
