/*

The model files that ship with the repository, as examples.

MSP 16 Apr 2022 first definitions

@IMPROVE: move all parameters for simulation into this file? I.e. clean up translate.h?

*/
#ifndef MODELFILES_H
#define MODELFILES_H

using std::string;


string indir = "../data/";
string runtag = "run9mlde";


// the virial time of the 'present day'
double reference_time = 0.0;

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
string sph_cache_name_mw   = indir+"SLGridSph.mw."+runtag;
string model_file_mw       = indir+"ErkalMW.model";
string coef_file_mw        = indir+"simpleoutcoef.nfac.mw."+runtag;
string orient_file_mw      = indir+"mw.orient."+runtag+".smth";

// The MW disc tables
// string cyl_cache_name_mw   = indir+"disc.cache."+runtag; // this is the maximum accuracy bais (too big for GitHub!)
string cyl_cache_name_mw   = indir+"disc.cache.compact."+runtag;
string cyl_coef_name_mw    = indir+"outcoef.disc."+runtag;
string cyl_orient_name_mw  = indir+"disc.orient."+runtag+".smth";

// The LMC halo tables
string sph_cache_name_lmc  = indir+"SLGridSph.lmc."+runtag;
string model_file_lmc      = indir+"ErkalLMC.model";
string coef_file_lmc       = indir+"simpleoutcoef.nfac.lmc."+runtag;
string orient_file_lmc     = indir+"lmc.orient."+runtag+".smth";

// inform the compiler that we have a model
#define MODELDEFINED 1

#endif
