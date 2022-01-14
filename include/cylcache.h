/*
sphcache.h

functions to handle the preparations for the cachefile(s)

clean version, MSP 5 May 2020
MSP 21 Dec 2021 update yaml headers


Note that cachetable reading is slow.
We might be able to find a better way using memmaps, which would also
make querying the tables faster?

Check on mapping values: are these okay?

 */
#ifndef CYLCACHE_H
#define CYLCACHE_H


// should we have some way to block this out just in case?
#if HAVEYAML
#include "yaml-cpp/yaml.h"	// YAML support
#endif

//using namespace std;
using std::cout, std::cerr, std::endl, std::setw;

// create 2-,3- and 4-d array types
typedef boost::multi_array<double, 4> array_type4;
typedef boost::multi_array<double, 3> array_type3;
typedef boost::multi_array<double, 2> array_type2;


// create the cylindrical cachefile array
struct CylCache
{
  // from the cachefile itself
  int MMAX;            // the number of azimuthal harmonics
  int NUMX;            // points in the radial direction
  int NUMY;            // points in the vertical direction
  int NMAX;            // maximum radial order for basis construction
  int NORDER;          // number of radial terms
  bool DENS;            // density flag
  bool CMAP;            // mapping flag
  bool CMAPR;            // mapping flag R
  bool CMAPZ;            // mapping flag Z
  double RMIN;         // the minimum expansion radius
  double RMAX;         // the maximum expansion radius
  double ASCALE;       // the scaling value for scale length
  double HSCALE;       // the scaling value for scale height
  double CYLMASS;      // the mass of the initial component
  double TNOW;         // the time of construction

  array_type4 potC;    // the cosine potential table,      sized MMAX,NORDER,NUMX+1,NUMY+1
  array_type4 potS;    // the   sine potential table,      sized MMAX,NORDER,NUMX+1,NUMY+1
  array_type4 rforceC; // the cosine radial force table,   sized MMAX,NORDER,NUMX+1,NUMY+1
  array_type4 rforceS; // the   sine radial table,         sized MMAX,NORDER,NUMX+1,NUMY+1
  array_type4 zforceC; // the cosine vertical force table, sized MMAX,NORDER,NUMX+1,NUMY+1
  array_type4 zforceS; // the   sine vertical force table, sized MMAX,NORDER,NUMX+1,NUMY+1
  array_type4 densC;   // the cosine density table,        sized MMAX,NORDER,NUMX+1,NUMY+1
  array_type4 densS;   // the   sine density table,        sized MMAX,NORDER,NUMX+1,NUMY+1


  double Rtable;        // maximum table radius

  double dX;           // the spacing in the scaled radius array
  double XMIN;         // the minimum scaled radius
  double XMAX;         // the maximum scaled radius

  double dY;           // the spacing in the scaled vertical array
  double YMIN;         // the minimum scaled vertical value
  double YMAX;         // the maximum scaled vertical value


};

struct CylForce
{

  array_type2 potC;
  array_type2 potS;

  array_type2 rforceC;
  array_type2 rforceS;

  array_type2 zforceC;
  array_type2 zforceS;

};


void read_cyl_cache(string& cyl_cache_name, CylCache& cachetable)
{

  ifstream in(cyl_cache_name.c_str());

  cerr << "cylcache.read_cyl_cache: trying to read cached table. . . ";


  // Attempt to read magic number
  //
  unsigned int tmagic;
  in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));

  //cout << setw(14) << tmagic << endl;

  if (tmagic == 202004385) {

#if FORMAT
    cout << "NEW FORMAT" << endl;
#endif

    unsigned ssize;
    in.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));

    // Make and read char buffer
    //
    auto buf = std::make_unique<char[]>(ssize+1);
    in.read(buf.get(), ssize);
    buf[ssize] = 0;		// Null terminate

#if HAVEYAML
    YAML::Node node;

    try {
      node = YAML::Load(buf.get());
    } catch (YAML::Exception& error) {
      std::cerr << "YAML: error parsing <" << buf.get() << "> "
		    << "in " << __FILE__ << ":" << __LINE__ << std::endl
		    << "YAML error: " << error.what() << std::endl;
	throw error;
    }

      // Get parameters
      //
      cachetable.MMAX    = node["mmax"  ].as<int>();
      cachetable.NUMX    = node["numx"  ].as<int>();
      cachetable.NUMY    = node["numy"  ].as<int>();
      cachetable.NMAX    = node["nmax"  ].as<int>();
      cachetable.NORDER  = node["norder"].as<int>();
      cachetable.DENS    = node["dens"  ].as<bool>();
      cachetable.RMIN    = node["rmin"  ].as<double>();
      cachetable.RMAX    = node["rmax"  ].as<double>();
      cachetable.ASCALE  = node["ascl"  ].as<double>();
      cachetable.HSCALE  = node["hscl"  ].as<double>();
      cachetable.CYLMASS = node["cmass" ].as<double>();
      cachetable.TNOW    = node["time"  ].as<double>();

      if (node["cmap"]) 	// Backwards compatibility
	      cachetable.CMAPR   = node["cmap" ].as<int>();
      else
	      cachetable.CMAPR   = node["cmapr"].as<int>();

      if (node["cmapz"])	// Backwards compatibility
	      cachetable.CMAPZ = node["cmapz"].as<int>();
      else
	      cachetable.CMAPZ = 1;//CMAPZ;

#else // compile without libyaml

  cachetable.CMAPZ = 1;

  std::string yamlblob = buf.get();

  // loop through string parsing
  // https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
  // but NOT the top vote getter!
  std::string delim = "\n";
  std::string entry;
  auto start = 0U;
  auto end = yamlblob.find(delim);
  std::string token2,token3;
  while (end != std::string::npos)
  {
      // retreive the entry
      entry = yamlblob.substr(start, end - start);

      // split the entry up to the ':'
      token2 = entry.substr(0, entry.find(":"));

      // split the entry after the ':'
      token3 = entry.substr(entry.find(":")+2,entry.length()-entry.find(":")-2);
      //std::cout<<token2<<"--"<<token3<<endl;

      if (token2.compare("mmax") == 0)   cachetable.MMAX    = std::stoi(token3);
      if (token2.compare("norder") == 0) cachetable.NORDER  = std::stoi(token3);
      if (token2.compare("numx") == 0)   cachetable.NUMX    = std::stoi(token3);
      if (token2.compare("numy") == 0)   cachetable.NUMY    = std::stoi(token3);
      if (token2.compare("nmax") == 0)   cachetable.NMAX    = std::stoi(token3);
      if (token2.compare("rmin") == 0)   cachetable.RMIN    = std::stod(token3);
      if (token2.compare("rmax") == 0)   cachetable.RMAX    = std::stod(token3);
      if (token2.compare("ascl") == 0)   cachetable.ASCALE  = std::stod(token3);
      if (token2.compare("hscl") == 0)   cachetable.HSCALE  = std::stod(token3);
      if (token2.compare("cmass") == 0)  cachetable.CYLMASS = std::stod(token3);
      if (token2.compare("time") == 0)   cachetable.TNOW    = std::stod(token3);
      if (token2.compare("dens") == 0)   cachetable.DENS    = token3.compare("true")==0;

      if (token2.compare("cmap") == 0)   cachetable.CMAPR   = std::stoi(token3);
      if (token2.compare("cmapr") == 0)  cachetable.CMAPR   = std::stoi(token3);
      if (token2.compare("cmapz") == 0)  cachetable.CMAPZ   = std::stoi(token3);

      // advance counters
      start = end + delim.length();
      end = yamlblob.find(delim, start);
  }

  // last one is always 'time' which we can ignore
  //std::cout << yamlblob.substr(start, end) << "|";


  //std::cout << cachetable.MMAX  << " " << cachetable.ASCALE  << std::endl;


#endif

  } else {

#if FORMAT
    cout << "OLD FORMAT" << endl;
#endif
    int tmp;

    // rewind to the beginning of the file
    //
    in.clear();
    in.seekg(0);

    // read the parameters from the cachefile
    in.read((char *)&cachetable.MMAX,   sizeof(int));
    in.read((char *)&cachetable.NUMX,   sizeof(int));
    in.read((char *)&cachetable.NUMY,   sizeof(int));
    in.read((char *)&cachetable.NMAX,   sizeof(int));
    in.read((char *)&cachetable.NORDER, sizeof(int));
    in.read((char *)&tmp,               sizeof(int)); if (tmp) cachetable.DENS=true;
    in.read((char *)&tmp,               sizeof(int)); if (tmp) cachetable.CMAP=true;
    in.read((char *)&cachetable.RMIN,   sizeof(double));
    in.read((char *)&cachetable.RMAX,   sizeof(double));
    in.read((char *)&cachetable.ASCALE, sizeof(double));
    in.read((char *)&cachetable.HSCALE, sizeof(double));
    in.read((char *)&cachetable.CYLMASS,sizeof(double));
    in.read((char *)&cachetable.TNOW,   sizeof(double));

  }

#if ORDERS
  cout << setw(14) << cachetable.MMAX << setw(14) << cachetable.NORDER << endl;
#endif

  // resize the arrays
  cachetable.potC.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER][cachetable.NUMX+1][cachetable.NUMY+1]);
  cachetable.potS.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER][cachetable.NUMX+1][cachetable.NUMY+1]);

  cachetable.rforceC.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER][cachetable.NUMX+1][cachetable.NUMY+1]);
  cachetable.rforceS.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER][cachetable.NUMX+1][cachetable.NUMY+1]);

  cachetable.zforceC.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER][cachetable.NUMX+1][cachetable.NUMY+1]);
  cachetable.zforceS.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER][cachetable.NUMX+1][cachetable.NUMY+1]);

  if (cachetable.DENS) {
    cachetable.densC.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER][cachetable.NUMX+1][cachetable.NUMY+1]);
    cachetable.densS.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER][cachetable.NUMX+1][cachetable.NUMY+1]);
  }


  // read in cosine components
  for (int m=0; m<=cachetable.MMAX; m++) {

    // read in eigenvector values
    for (int n=0; n<cachetable.NORDER; n++) {

      for (int j=0; j<=cachetable.NUMX; j++) {
	for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.potC[m][n][j][k], sizeof(double));
      }

      for (int j=0; j<=cachetable.NUMX; j++) {
	for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.rforceC[m][n][j][k], sizeof(double));
      }

      for (int j=0; j<=cachetable.NUMX; j++) {
	for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.zforceC[m][n][j][k], sizeof(double));
      }

      if (cachetable.DENS) {
        for (int j=0; j<=cachetable.NUMX; j++) {
	  for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.densC[m][n][j][k], sizeof(double));
	}
      }
    } // end NORDER loop
  } // end MMAX loop

  // read in sine components (note: no 0 order for sine terms)
  for (int m=1; m<=cachetable.MMAX; m++) {

    // read in eigenvector values
    for (int n=0; n<cachetable.NORDER; n++) {

      for (int j=0; j<=cachetable.NUMX; j++) {
	for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.potS[m][n][j][k], sizeof(double));
      }

      for (int j=0; j<=cachetable.NUMX; j++) {
	for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.rforceS[m][n][j][k], sizeof(double));
      }

      for (int j=0; j<=cachetable.NUMX; j++) {
	for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.zforceS[m][n][j][k], sizeof(double));
      }

      if (cachetable.DENS) {
        for (int j=0; j<=cachetable.NUMX; j++) {
	  for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.densS[m][n][j][k], sizeof(double));
	}
      }
    } // end NORDER loop
  } // end MMAX loop


  // set the table parameters
    //set_table_params
    //    calculate scaled boundary values for the parameter table

  cachetable.Rtable  = M_SQRT1_2 * cachetable.RMAX;

  // calculate radial scalings
  cachetable.XMIN    = r_to_xi(cachetable.RMIN*cachetable.ASCALE,cachetable.CMAP,cachetable.ASCALE);
  cachetable.XMAX    = r_to_xi(cachetable.Rtable*cachetable.ASCALE,cachetable.CMAP,cachetable.ASCALE);
  cachetable.dX      = (cachetable.XMAX - cachetable.XMIN)/cachetable.NUMX;

  // calculate vertical scalings
  cachetable.YMIN    = z_to_y(-cachetable.Rtable*cachetable.ASCALE,cachetable.HSCALE);
  cachetable.YMAX    = z_to_y( cachetable.Rtable*cachetable.ASCALE,cachetable.HSCALE);
  cachetable.dY      = (cachetable.YMAX - cachetable.YMIN)/cachetable.NUMY;


  std::cerr << "success!!" << std::endl;
}



void get_pot(double r, double z, CylCache cachetable, array_type2& Vc, array_type2& Vs)
{
  Vc.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER]);
  Vs.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER]);

  if (z/cachetable.ASCALE > cachetable.Rtable) z =  cachetable.Rtable*cachetable.ASCALE;
  if (z/cachetable.ASCALE <-cachetable.Rtable) z = -cachetable.Rtable*cachetable.ASCALE;

  double X = (r_to_xi(r,cachetable.CMAP,cachetable.ASCALE) - cachetable.XMIN)/cachetable.dX;
  double Y = (z_to_y(z,cachetable.HSCALE) - cachetable.YMIN)/cachetable.dY;

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

      Vc[mm][n] =
	(
	 cachetable.potC[mm][n][ix  ][iy  ] * c00 +
	 cachetable.potC[mm][n][ix+1][iy  ] * c10 +
	 cachetable.potC[mm][n][ix  ][iy+1] * c01 +
	 cachetable.potC[mm][n][ix+1][iy+1] * c11
	 );

      // get sine values for m>0
      if (mm) {

	Vs[mm][n] =
	  (
	   cachetable.potS[mm][n][ix  ][iy  ] * c00 +
	   cachetable.potS[mm][n][ix+1][iy  ] * c10 +
	   cachetable.potS[mm][n][ix  ][iy+1] * c01 +
	   cachetable.potS[mm][n][ix+1][iy+1] * c11
	   );
      }

    } // end NORDER loop
  } // end MMAX loop

}



void get_table_forces(double r, double z, CylCache cachetable, CylForce& forcetable)
{

  // return 2d tables required to compute the forces

  forcetable.potC.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER]);
  forcetable.potS.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER]);

  forcetable.rforceC.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER]);
  forcetable.rforceS.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER]);

  forcetable.zforceC.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER]);
  forcetable.zforceS.resize(boost::extents[cachetable.MMAX+1][cachetable.NORDER]);


  if (z/cachetable.ASCALE > cachetable.Rtable) z =  cachetable.Rtable*cachetable.ASCALE;
  if (z/cachetable.ASCALE <-cachetable.Rtable) z = -cachetable.Rtable*cachetable.ASCALE;

  double X = (r_to_xi(r,cachetable.CMAP,cachetable.ASCALE) - cachetable.XMIN)/cachetable.dX;
  double Y = (z_to_y(z,cachetable.HSCALE) - cachetable.YMIN)/cachetable.dY;

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

      forcetable.potC[mm][n] =
	(
	 cachetable.potC[mm][n][ix  ][iy  ] * c00 +
	 cachetable.potC[mm][n][ix+1][iy  ] * c10 +
	 cachetable.potC[mm][n][ix  ][iy+1] * c01 +
	 cachetable.potC[mm][n][ix+1][iy+1] * c11
	 );

      forcetable.rforceC[mm][n] =
	(
	 cachetable.rforceC[mm][n][ix  ][iy  ] * c00 +
	 cachetable.rforceC[mm][n][ix+1][iy  ] * c10 +
	 cachetable.rforceC[mm][n][ix  ][iy+1] * c01 +
	 cachetable.rforceC[mm][n][ix+1][iy+1] * c11
	 );

      forcetable.zforceC[mm][n] =
	(
	 cachetable.zforceC[mm][n][ix  ][iy  ] * c00 +
	 cachetable.zforceC[mm][n][ix+1][iy  ] * c10 +
	 cachetable.zforceC[mm][n][ix  ][iy+1] * c01 +
	 cachetable.zforceC[mm][n][ix+1][iy+1] * c11
	 );

      // get sine values for m>0
      if (mm) {

	forcetable.potS[mm][n] =
	  (
	   cachetable.potS[mm][n][ix  ][iy  ] * c00 +
	   cachetable.potS[mm][n][ix+1][iy  ] * c10 +
	   cachetable.potS[mm][n][ix  ][iy+1] * c01 +
	   cachetable.potS[mm][n][ix+1][iy+1] * c11
	   );

        forcetable.rforceS[mm][n] =
	  (
	   cachetable.rforceS[mm][n][ix  ][iy  ] * c00 +
	   cachetable.rforceS[mm][n][ix+1][iy  ] * c10 +
	   cachetable.rforceS[mm][n][ix  ][iy+1] * c01 +
	   cachetable.rforceS[mm][n][ix+1][iy+1] * c11
	   );

        forcetable.zforceS[mm][n] =
	  (
	   cachetable.zforceS[mm][n][ix  ][iy  ] * c00 +
	   cachetable.zforceS[mm][n][ix+1][iy  ] * c10 +
	   cachetable.zforceS[mm][n][ix  ][iy+1] * c01 +
	   cachetable.zforceS[mm][n][ix+1][iy+1] * c11
	   );
      }

    } // end NORDER loop
  } // end MMAX loop

}




/*
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

  // this step would be great to change, for speed purposes.
  double pval = cachetable.modeltable.pspline(r);

  for (int l=0; l<=cachetable.LMAX; l++) {
    for (int n=0; n<cachetable.NMAX; n++) {
      pottable[l][n] = (x1*cachetable.eftable[l][n][indx] + x2*cachetable.eftable[l][n][indx+1])/
      sqrt(cachetable.evtable[l][n]) * pval;
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


void get_dpotl(double r, SphCache& cachetable, array_type2& potd, array_type2& dpot) {

  // r comes in as the actual radius, NOT xi
  get_pot  (r, cachetable, potd);
  get_force(r, cachetable, dpot);

}

*/


#endif
