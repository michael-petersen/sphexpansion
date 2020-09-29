/*
sphcache.h

functions to handle the preparations for the cachefile(s)

MSP 22 Apr 2020 clean version
MSP 26 Sep 2020 add density calls: to be fully verified
MSP 27 Sep 2020 add USE_TABLE preprocessor directive


 */

// include the modelfile stuff
#include "sphmodel.h"

using namespace std;

// create 2- and 3-d array types
typedef boost::multi_array<double, 3> array_type3;
typedef boost::multi_array<double, 2> array_type2;

typedef vector<double> dbl_vector;
typedef boost::multi_array<vector<dbl_vector>, 2> efarray;

bool debug = false;

// use tabulated values for density rather than spline
#define USE_TABLE 1



// create the spherical cachefile array
struct SphCache
{
  // from the cachefile itself
  int LMAX;            // the number of azimuthal harmonics
  int NMAX;            // the number of radial terms
  int NUMR;            // the number of points per radial term
  int CMAP;            // mapping flag
  double RMIN;         // the minimum expansion radius
  double RMAX;         // the maximum expansion radius
  double SCL;          // the scaling value for radius
  array_type3 eftable; // the eigenfunction table, sized LMAX+1,NMAX,NUMR
  array_type2 evtable; // the eigenvector table, sized LMAX+1,NMAX
  efarray eftable2;    // test eftable version

  // added after the SphModel is read in, in init_table
  vector<double> r;    // the radius array, len NUMR
  vector<double> d0;   // the density array, len NUMR
  vector<double> p0;   // the potential array, len NUMR
  vector<double> xi;   // the scaled radius array, len NUMR
  double dxi;          // the spacing in the scaled radius array
  double xmin;         // the minimum scaled radius
  double xmax;         // the maximum scaled radius

  // make the whole table available, whatever
  SphModel modeltable;
  
};

void read_sph_cache (string& sph_cache_name, SphCache& cachetable) {

  ifstream in(sph_cache_name.c_str());
  if (!in) {
        cout << "sphcache::read_sph_cache: Unable to open file!\n";
        exit(1);
  }

  cerr << "sphcache.read_sph_cache: trying to read cached table . . . ";

  // read the parameters from the cachefile
  in.read((char *)&cachetable.LMAX, sizeof(int));
  in.read((char *)&cachetable.NMAX, sizeof(int));
  in.read((char *)&cachetable.NUMR, sizeof(int));
  in.read((char *)&cachetable.CMAP, sizeof(int));	
  in.read((char *)&cachetable.RMIN, sizeof(double));
  in.read((char *)&cachetable.RMAX, sizeof(double));
  in.read((char *)&cachetable.SCL,  sizeof(double));

  // debug
  if (debug) 
  cout << setw(10) << cachetable.LMAX
       << setw(10) << cachetable.NMAX
       << setw(10) << cachetable.NUMR
       << setw(10) << cachetable.CMAP
       << setw(10) << cachetable.RMIN
       << setw(10) << cachetable.RMAX
       << setw(10) << cachetable.SCL
       << endl;

  // resize the arrays
  cachetable.eftable.resize(boost::extents[cachetable.LMAX+1][cachetable.NMAX][cachetable.NUMR]);
  cachetable.evtable.resize(boost::extents[cachetable.LMAX+1][cachetable.NMAX]);
  
  cachetable.eftable2.resize(boost::extents[cachetable.LMAX+1][cachetable.NMAX]);

  int dummy;
  for (int l=0; l<=cachetable.LMAX; l++) {

  // read in the table element first: not necessary to save, vestigial from exp
  in.read((char *)&dummy, sizeof(int));

  // read in eigenvector values
  for (int n=0; n<cachetable.NMAX; n++) in.read((char *)&cachetable.evtable[l][n], sizeof(double));

  // read in eigenfunction values
  for (int n=0; n<cachetable.NMAX; n++) {
    for (int i=0; i<cachetable.NUMR; i++)
      in.read((char *)&cachetable.eftable[l][n][i], sizeof(double));
    }
  }

  //cout << cachetable.RMIN << endl;
  std::cerr << "success!!" << std::endl;
}



void init_table(SphModel& sphmodel, SphCache& cachetable)
{
  // duplicate the whole sphmodel for access:
  cachetable.modeltable = sphmodel;
  // is this too painful for the memory footprint?
  
  // add the relevant vectors to the cachetable from the sphmodel
  cachetable.xi.resize(cachetable.NUMR);
  cachetable.r.resize(cachetable.NUMR);
  cachetable.p0.resize(cachetable.NUMR);
  cachetable.d0.resize(cachetable.NUMR);

  if (cachetable.CMAP==1) {
    cachetable.xmin = (cachetable.RMIN/cachetable.SCL - 1.0)/(cachetable.RMIN/cachetable.SCL + 1.0);
    cachetable.xmax = (cachetable.RMAX/cachetable.SCL - 1.0)/(cachetable.RMAX/cachetable.SCL + 1.0);
  }
  else if (cachetable.CMAP==2) {
    cachetable.xmin = log(cachetable.RMIN);
    cachetable.xmax = log(cachetable.RMAX);
  } else {
    cachetable.xmin = cachetable.RMIN;
    cachetable.xmax = cachetable.RMAX;
  }
  cachetable.dxi = (cachetable.xmax-cachetable.xmin)/(cachetable.NUMR-1);

  // rebuild to the same grid as the eigenfunctions
  for (int i=0; i<cachetable.NUMR; i++) {
    cachetable.xi[i] = cachetable.xmin + cachetable.dxi*i;
    cachetable.r[i]  = xi_to_r(cachetable.xi[i],cachetable.CMAP,cachetable.SCL);
    cachetable.p0[i] = get_pot(sphmodel, cachetable.r[i]);
    cachetable.d0[i] = get_density(sphmodel, cachetable.r[i]);
  }

}




void get_pot(double& r, SphCache& cachetable, array_type2& pottable)
{
  pottable.resize(boost::extents[cachetable.LMAX+1][cachetable.NMAX]);

  double xi;
  xi = r_to_xi(r, cachetable.CMAP, cachetable.SCL);
    
  if (cachetable.CMAP==1) {
        if (xi<-1.0) xi=-1.0;
        if (xi>=1.0) xi=1.0-1.0e-08;
  }

  int indx = (int)( (xi-cachetable.xmin)/cachetable.dxi );
  if (indx<0) indx = 0;
  if (indx>cachetable.NUMR-2) indx = cachetable.NUMR - 2;

  double x1 = (cachetable.xi[indx+1] - xi)/cachetable.dxi;
  double x2 = (xi - cachetable.xi[indx])/cachetable.dxi;

  for (int l=0; l<=cachetable.LMAX; l++) {
    for (int n=0; n<cachetable.NMAX; n++) {

      #ifdef USETABLE      
      pottable[l][n] = (x1*cachetable.eftable[l][n][indx] + x2*cachetable.eftable[l][n][indx+1])/
	sqrt(cachetable.evtable[l][n]) * (x1*cachetable.p0[indx] + x2*cachetable.p0[indx+1]);
      #else
      pottable[l][n] = (x1*cachetable.eftable[l][n][indx] + x2*cachetable.eftable[l][n][indx+1])/
	sqrt(cachetable.evtable[l][n]) * cachetable.modeltable.pspline(r);
      #endif
    }
  }
}


void get_force(double& r, SphCache& cachetable, array_type2& forcetable) {

  // see the equivalent call, get_force in SLGridMP2.cc
  
  // must have already run init_table

  forcetable.resize(boost::extents[cachetable.LMAX+1][cachetable.NMAX]);

  double xi;
  xi = r_to_xi(r, cachetable.CMAP, cachetable.SCL);
    
  if (cachetable.CMAP==1) {
        if (xi<-1.0) xi=-1.0;
        if (xi>=1.0) xi=1.0-1.0e-08;
  }

  int indx = (int)( (xi-cachetable.xmin)/cachetable.dxi );
  if (indx<1) indx = 1;
  if (indx>cachetable.NUMR-2) indx = cachetable.NUMR - 2;

  double p = (xi - cachetable.xi[indx])/cachetable.dxi;
  
				// Use three point formula

				// Point -1: indx-1
				// Point  0: indx
				// Point  1: indx+1

  for (int l=0; l<=cachetable.LMAX; l++) {
    for (int n=0; n<cachetable.NMAX; n++) {
      forcetable[l][n] = d_xi_to_r(xi,cachetable.CMAP,cachetable.SCL)/cachetable.dxi * (
			     (p - 0.5)*cachetable.eftable[l][n][indx-1]*cachetable.p0[indx-1]
			     -2.0*p*cachetable.eftable[l][n][indx]*cachetable.p0[indx]
			     + (p + 0.5)*cachetable.eftable[l][n][indx+1]*cachetable.p0[indx+1]
			     ) / sqrt(cachetable.evtable[l][n]);
    }
  }

}

void get_density(double& r, SphCache& cachetable, array_type2& densitytable)
{

  // see equivalent call in SLGridMP2.cc
  
  densitytable.resize(boost::extents[cachetable.LMAX+1][cachetable.NMAX]);

  double xi;
  xi = r_to_xi(r, cachetable.CMAP, cachetable.SCL);
    
  if (cachetable.CMAP==1) {
        if (xi<-1.0) xi=-1.0;
        if (xi>=1.0) xi=1.0-1.0e-08;
  }

  int indx = (int)( (xi-cachetable.xmin)/cachetable.dxi );
  if (indx<0) indx = 0;
  if (indx>cachetable.NUMR-2) indx = cachetable.NUMR - 2;

  double x1 = (cachetable.xi[indx+1] - xi)/cachetable.dxi;
  double x2 = (xi - cachetable.xi[indx])/cachetable.dxi;

  //cout << cachetable.d0[indx] << endl;


  for (int l=0; l<=cachetable.LMAX; l++) {
    for (int n=0; n<cachetable.NMAX; n++) {

      #ifdef USETABLE
      densitytable[l][n] = (x1*cachetable.eftable[l][n][indx] + x2*cachetable.eftable[l][n][indx+1]) *
        sqrt(cachetable.evtable[l][n]) * (x1*cachetable.d0[indx] + x2*cachetable.d0[indx+1]);
      #else
      densitytable[l][n] = (x1*cachetable.eftable[l][n][indx] + x2*cachetable.eftable[l][n][indx+1]) *
        sqrt(cachetable.evtable[l][n]) * cachetable.modeltable.dspline(r);
      #endif
    }
  }
}


void get_dpotl(double r, SphCache& cachetable, array_type2& potd, array_type2& dpot)
{
  // r comes in as the actual radius, NOT xi
  get_pot  (r, cachetable, potd);
  get_force(r, cachetable, dpot);
}

void get_dpotl_density(double r, SphCache& cachetable, array_type2& potd, array_type2& dpot, array_type2& dend)
{
  // r comes in as the actual radius, NOT xi
  get_pot  (r, cachetable, potd);
  get_force(r, cachetable, dpot);
  get_density(r, cachetable, dend);
}
