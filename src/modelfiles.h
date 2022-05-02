/*

The model files that ship with the repository, as examples.

MSP 16 Apr 2022 first definitions

*/
#ifndef MODELFILES_H
#define MODELFILES_H

using std::string;

std::string get_datadir (const std::string& str)
{
  std::size_t found = str.find_last_of("/\\");
  return (str.substr(0,found+1)).append("../data/");
}

string model_datadir = get_datadir(__FILE__);
string runtag = "run9mlde";

// the virial time of the 'present day'
double reference_time = 1.19925; // (in virial)
double native_timestep = 0.000125; // native timestep of the simulation (in virial)
int    native_steps = reference_time/native_timestep; // number of native timesteps in the simulation (does this match NUMT?)

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

// The MW halo tables
string sph_cache_name_mw   = model_datadir+"SLGridSph.mw."+runtag;
string model_file_mw       = model_datadir+"ErkalMW.model";
string coef_file_mw        = model_datadir+"simpleoutcoef.nfac.mw."+runtag;
string orient_file_mw      = model_datadir+"mw.orient."+runtag+".smth";

// The MW disc tables
string cyl_cache_name_mw   = model_datadir+"disc.cache.compact."+runtag;
string cyl_coef_name_mw    = model_datadir+"outcoef.disc."+runtag;
string cyl_orient_name_mw  = model_datadir+"disc.orient."+runtag+".smth";

// The LMC halo tables
string sph_cache_name_lmc  = model_datadir+"SLGridSph.lmc."+runtag;
string model_file_lmc      = model_datadir+"ErkalLMC.model";
string coef_file_lmc       = model_datadir+"simpleoutcoef.nfac.lmc."+runtag;
string orient_file_lmc     = model_datadir+"lmc.orient."+runtag+".smth";

// inform the compiler that we have a model
#define MODELDEFINED 1

#endif
