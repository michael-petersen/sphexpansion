/*
translate.h

functions to take virial units and make physical units

notes from call
-1 Gyr ~ 1 kpc/km/s ... but we'll stick with kpc/km/s

MSP 23 Apr 2020 clean version
MSP  6 May 2020 close the loop with Msun addition
MSP 28 Sep 2020 add density translation
MSP 14 May 2021 added header guard
MSP 16 Apr 2022 move simulation-specific definitions to different header
MSP 19 Apr 2022 add potential translation

@IMPROVE: write overloaded versions of the functions

 */
#ifndef TRANSLATE_H
#define TRANSLATE_H
using std::cout, std::cerr, std::endl, std::setw;

/*
for translations to work, the following quantities must be defined:
mw_vel_scale
mw_time_scale
mw_force_scale
mw_mass_scale

if not defined, all are assumed to be equal to unity (virial units)
*/

#if MODELDEFINED == 0
double mw_vel_scale            = 1.
double mw_time_scale           = 1.
double mw_force_scale          = 1.
double mw_mass_scale           = 1.
#endif


void virial_to_physical_density(double densvir, double& densphys)
{
  // convert Mvir/Rvir^3 to Msun/pc^3

  densphys = densvir * mw_mass_scale / pow(1000*mw_virial_radius,3);
}

void virial_to_physical_time(double tvir, double& tphys)
{
  tphys = mw_time_scale * tvir;
}

double virial_to_physical_time_return(double tvir)
{
  return mw_time_scale * tvir;
}

void virial_to_physical_potential(double pvir, double& pphys)
{
  pphys = mw_mass_scale * mw_force_scale * pvir;
}

void virial_to_physical_force(double  fxvir, double   fyvir, double   fzvir,
			      double& fxphys, double& fyphys, double& fzphys)
{
  fxphys = mw_force_scale * fxvir;
  fyphys = mw_force_scale * fyvir;
  fzphys = mw_force_scale * fzvir;
}

double virial_to_physical_force(double  fxvir)
{
  return mw_force_scale * fxvir;
}

void virial_to_physical_length(double xvir, double   yvir, double   zvir,
			      double& xphys, double& yphys, double& zphys)
{
  xphys = mw_virial_radius * xvir;
  yphys = mw_virial_radius * yvir;
  zphys = mw_virial_radius * zvir;
}

double virial_to_physical_length(double xvir)
{
  // no overload, single coordinate transformation
  return mw_virial_radius * xvir;
}

void virial_to_physical_velocity(double vxvir, double vyvir, double vzvir,
			double& vxphys, double& vyphys, double& vzphys)
{
  vxphys = mw_vel_scale * vxvir;
  vyphys = mw_vel_scale * vyvir;
  vzphys = mw_vel_scale * vzvir;
}

double virial_to_physical_velocity(double vxvir)
{
  return mw_vel_scale * vxvir;
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

double physical_to_virial_time_return(double tphys)
{
  return tphys/mw_time_scale;
}

void physical_to_virial_potential(double pphys, double& pvir)
{
  pvir  = pphys / mw_mass_scale / mw_force_scale;
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

double physical_to_virial_length(double xphys)
{
  return xphys/mw_virial_radius;
}

void physical_to_virial_velocity(double  vxphys, double vyphys, double vzphys,
		         	 double& vxvir, double& vyvir, double& vzvir)
{
  vxvir = vxphys/mw_vel_scale;
  vyvir = vyphys/mw_vel_scale;
  vzvir = vzphys/mw_vel_scale;
}

double physical_to_virial_velocity(double  vxphys)
{
  return vxphys/mw_vel_scale;
}

#endif
