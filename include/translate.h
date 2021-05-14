/*
translate.h

functions to take virial units and make physical units

notes from call
-1 Gyr ~ 1 kpc/km/s ... but we'll stick with kpc/km/s

MSP 23 Apr 2020 clean version
MSP  6 May 2020 close the loop with Msun addition
MSP 28 Sep 2020 add density translation
MSP 14 May 2021 added header guard


 */
#ifndef TRANSLATE_H
#define TRANSLATE_H

using namespace std;

// set these parameters to tune the simulation virial units to the MW
double mw_virial_radius        = 282.;                                  // virial radius of the MW, kpc/Rvir
double solar_circular_velocity = 190.;                                  // circular velocity at the solar circle, km/s
double rotation_peak           = 1.4 ;                                  // peak of the rotation curve from the MW spherical model (Vvir)

// derived quantities
double astronomicalG           = 0.0000043009125;                       // gravitational constant, (km/s)^2 * kpc / Msun
double mw_vel_scale            = solar_circular_velocity/rotation_peak; // velocity scale of the MW, km/s/Vvir
double mw_time_scale           = mw_virial_radius/mw_vel_scale;         // time scale for the MW, kpc/km/s
double mw_force_scale          = mw_vel_scale/mw_time_scale;            // in km/s/s
double mw_mass_scale           = mw_vel_scale * mw_vel_scale *
                                 mw_virial_radius / astronomicalG;       // mass units of the simulation, Msun


void virial_to_physical_density(double densvir, double& densphys)
{
  // convert Mvir/Rvir^3 to Msun/pc^3
    
  densphys = densvir * mw_mass_scale / pow(1000*mw_virial_radius,3);
}

void virial_to_physical_time(double tvir, double& tphys)
{
  tphys = mw_time_scale * tvir;
}

void virial_to_physical_force(double  fxvir, double   fyvir, double   fzvir,
			      double& fxphys, double& fyphys, double& fzphys)
{
  fxphys = mw_force_scale * fxvir;
  fyphys = mw_force_scale * fyvir;
  fzphys = mw_force_scale * fzvir;
}

void virial_to_physical_length(double xvir, double   yvir, double   zvir,
			      double& xphys, double& yphys, double& zphys)
{
  xphys = mw_virial_radius * xvir;
  yphys = mw_virial_radius * yvir;
  zphys = mw_virial_radius * zvir;
}

void virial_to_physical_velocity(double vxvir, double vyvir, double vzvir,
			double& vxphys, double& vyphys, double& vzphys)
{
  vxphys = mw_vel_scale * vxvir;
  vyphys = mw_vel_scale * vyvir;
  vzphys = mw_vel_scale * vzvir;
}

void physical_to_virial_density(double densphys, double& densvir)
{
  // convert Msun/pc^3 to Mvir/Rvir^3 
    
  densvir = densphys * pow(1000*mw_virial_radius,3)/mw_mass_scale;
}

void physical_to_virial_time(double tphys, double& tvir)
{
  tvir = tphys/mw_time_scale;
}

void physical_to_virial_force(double  fxphys, double fyphys, double fzphys,
			      double& fxvir, double& fyvir, double& fzvir)
{
  fxvir = fxphys/mw_force_scale;
  fyvir = fyphys/mw_force_scale;
  fzvir = fzphys/mw_force_scale;
}

void physical_to_virial_length(double xphys, double yphys, double zphys,
			      double& xvir, double& yvir, double& zvir)
{
  xvir = xphys/mw_virial_radius;
  yvir = yphys/mw_virial_radius;
  zvir = zphys/mw_virial_radius;
}

void physical_to_virial_velocity(double  vxphys, double vyphys, double vzphys,
		         	 double& vxvir, double& vyvir, double& vzvir)
{
  vxvir = vxphys/mw_vel_scale;
  vyvir = vyphys/mw_vel_scale;
  vzvir = vzphys/mw_vel_scale;
}


#endif
