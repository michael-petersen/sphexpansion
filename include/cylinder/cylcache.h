/*
cylcache.h

functions to handle the preparations for the cylindrical cachefile(s) from EXP

MSP  5 May 2020 clean version
MSP 21 Dec 2021 update yaml headers
MSP 18 Apr 2022 update CMAPR,CMAPZ parameters

Note that cachetable reading is slow.
We might be able to find a better way using memmaps, which would also make querying the tables faster?

 */
#ifndef CYLCACHE_H
#define CYLCACHE_H



using std::cout, std::cerr, std::endl, std::setw, std::vector, std::ifstream, std::ios, std::string, std::ofstream, std::istringstream;

// Eigen MatrixXd, std::vector <MatrixXd>
#include <Eigen/StdVector>
#include <Eigen/Dense>
using Eigen::MatrixXd;

typedef std::vector< std::vector< MatrixXd > > array_type4;


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
  int CMAPZ;            // mapping flag Z
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

  MatrixXd potC;
  MatrixXd potS;

  MatrixXd rforceC;
  MatrixXd rforceS;

  MatrixXd zforceC;
  MatrixXd zforceS;

};


void read_cyl_cache(string& cyl_cache_name, CylCache& cachetable)
{

  ifstream in(cyl_cache_name.c_str());

  if (!in) {
    throw "cylcache::read_cyl_cache: Unable to open file!\n";
  }

  std::cerr << "cylcache.read_cyl_cache: trying to read cached table. . . ";


  // Attempt to read magic number
  //
  unsigned int tmagic;
  in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));
  //in.read((char *)(&tmagic), sizeof(unsigned int));

#if DEBUGCACHE
  cout << "FORMAT FLAG=" << setw(14) << tmagic << " (vs 202004385) at " << in.tellg() << " for size " << sizeof(unsigned int) << endl;
#endif

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
    buf[ssize] = 0;    // Null terminate


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

    cachetable.CMAP = cachetable.CMAPR;

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

    cachetable.CMAPR = cachetable.CMAP;
    cachetable.CMAPZ = 1; // always 1 in old style

  }

#if ORDERS
  cout << setw(14) << cachetable.MMAX << setw(14) << cachetable.NORDER << endl;
  cout << setw(14) << cachetable.CMAPR << setw(14) << cachetable.CMAP << setw(14) << cachetable.CMAPZ << endl;
#endif

  // resize the arrays
  cachetable.potC.resize(cachetable.MMAX+1);
  cachetable.potS.resize(cachetable.MMAX+1);

  cachetable.rforceC.resize(cachetable.MMAX+1);
  cachetable.rforceS.resize(cachetable.MMAX+1);

  cachetable.zforceC.resize(cachetable.MMAX+1);
  cachetable.zforceS.resize(cachetable.MMAX+1);

  if (cachetable.DENS) {
    cachetable.densC.resize(cachetable.MMAX+1);
    cachetable.densS.resize(cachetable.MMAX+1);
  }


  // read in cosine components
  for (int m=0; m<=cachetable.MMAX; m++) {

    cachetable.potC[m].resize(cachetable.NORDER);
    cachetable.potS[m].resize(cachetable.NORDER);

    cachetable.rforceC[m].resize(cachetable.NORDER);
    cachetable.rforceS[m].resize(cachetable.NORDER);

    cachetable.zforceC[m].resize(cachetable.NORDER);
    cachetable.zforceS[m].resize(cachetable.NORDER);

    if (cachetable.DENS) {
      cachetable.densC[m].resize(cachetable.NORDER);
      cachetable.densS[m].resize(cachetable.NORDER);
    }

    // read in eigenvector values
    for (int n=0; n<cachetable.NORDER; n++) {

      cachetable.potC[m][n].resize(cachetable.NUMX+1,cachetable.NUMY+1);
      cachetable.potS[m][n].resize(cachetable.NUMX+1,cachetable.NUMY+1);

      cachetable.rforceC[m][n].resize(cachetable.NUMX+1,cachetable.NUMY+1);
      cachetable.rforceS[m][n].resize(cachetable.NUMX+1,cachetable.NUMY+1);

      cachetable.zforceC[m][n].resize(cachetable.NUMX+1,cachetable.NUMY+1);
      cachetable.zforceS[m][n].resize(cachetable.NUMX+1,cachetable.NUMY+1);

      if (cachetable.DENS) {
        cachetable.densC[m][n].resize(cachetable.NUMX+1,cachetable.NUMY+1);
        cachetable.densS[m][n].resize(cachetable.NUMX+1,cachetable.NUMY+1);
      }

      for (int j=0; j<=cachetable.NUMX; j++) {
        for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.potC[m][n](j,k), sizeof(double));
      }

      for (int j=0; j<=cachetable.NUMX; j++) {
        for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.rforceC[m][n](j,k), sizeof(double));
      }

      for (int j=0; j<=cachetable.NUMX; j++) {
         for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.zforceC[m][n](j,k), sizeof(double));
      }

      if (cachetable.DENS) {
        for (int j=0; j<=cachetable.NUMX; j++) {
           for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.densC[m][n](j,k), sizeof(double));
        }
      } // end if DENS
    } // end NORDER loop
  } // end MMAX loop

  // read in sine components (note: no 0 order for sine terms)
  for (int m=1; m<=cachetable.MMAX; m++) {

    // read in eigenvector values
    for (int n=0; n<cachetable.NORDER; n++) {

      for (int j=0; j<=cachetable.NUMX; j++) {
         for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.potS[m][n](j,k), sizeof(double));
      }

      for (int j=0; j<=cachetable.NUMX; j++) {
         for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.rforceS[m][n](j,k), sizeof(double));
      }

      for (int j=0; j<=cachetable.NUMX; j++) {
         for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.zforceS[m][n](j,k), sizeof(double));
      }

      if (cachetable.DENS) {
        for (int j=0; j<=cachetable.NUMX; j++) {
           for (int k=0; k<=cachetable.NUMY; k++) in.read((char *)&cachetable.densS[m][n](j,k), sizeof(double));
        }
      }
    } // end NORDER loop
  } // end MMAX loop


  // set the table parameters
    //set_table_params
    //    calculate scaled boundary values for the parameter table

  cachetable.Rtable  = M_SQRT1_2 * cachetable.RMAX;

  // calculate radial scalings
  cachetable.XMIN    = r_to_xi_cyl(cachetable.RMIN*cachetable.ASCALE,cachetable.CMAPR,cachetable.ASCALE);
  cachetable.XMAX    = r_to_xi_cyl(cachetable.Rtable*cachetable.ASCALE,cachetable.CMAPR,cachetable.ASCALE);
  cachetable.dX      = (cachetable.XMAX - cachetable.XMIN)/cachetable.NUMX;

  // calculate vertical scalings
  cachetable.YMIN    = z_to_y_cyl(-cachetable.Rtable*cachetable.ASCALE,cachetable.CMAPZ,cachetable.HSCALE);
  cachetable.YMAX    = z_to_y_cyl( cachetable.Rtable*cachetable.ASCALE,cachetable.CMAPZ,cachetable.HSCALE);
  cachetable.dY      = (cachetable.YMAX - cachetable.YMIN)/cachetable.NUMY;


  std::cerr << "success!!" << std::endl;
}



void get_pot(double r, double z, CylCache cachetable, MatrixXd& Vc, MatrixXd& Vs)
{
  Vc.resize(cachetable.MMAX+1,cachetable.NORDER);
  Vs.resize(cachetable.MMAX+1,cachetable.NORDER);

  if (z/cachetable.ASCALE > cachetable.Rtable) z =  cachetable.Rtable*cachetable.ASCALE;
  if (z/cachetable.ASCALE <-cachetable.Rtable) z = -cachetable.Rtable*cachetable.ASCALE;

  double X = (r_to_xi_cyl(r,cachetable.CMAPR,cachetable.ASCALE) - cachetable.XMIN)/cachetable.dX;
  double Y = (z_to_y_cyl(z,cachetable.CMAPZ,cachetable.HSCALE) - cachetable.YMIN)/cachetable.dY;

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

      Vc(mm,n) =
        (
         cachetable.potC[mm][n](ix  ,iy  ) * c00 +
         cachetable.potC[mm][n](ix+1,iy  ) * c10 +
         cachetable.potC[mm][n](ix  ,iy+1) * c01 +
         cachetable.potC[mm][n](ix+1,iy+1) * c11
         );

      // get sine values for m>0
      if (mm) {

        Vs(mm,n) =
          (
           cachetable.potS[mm][n](ix  ,iy  ) * c00 +
           cachetable.potS[mm][n](ix+1,iy  ) * c10 +
           cachetable.potS[mm][n](ix  ,iy+1) * c01 +
           cachetable.potS[mm][n](ix+1,iy+1) * c11
           );
            }

    } // end NORDER loop
  } // end MMAX loop

}



void get_table_forces(double r, double z, CylCache cachetable, CylForce& forcetable)
{

  // return 2d tables required to compute the forces

  forcetable.potC.resize(cachetable.MMAX+1,cachetable.NORDER);
  forcetable.potS.resize(cachetable.MMAX+1,cachetable.NORDER);

  forcetable.rforceC.resize(cachetable.MMAX+1,cachetable.NORDER);
  forcetable.rforceS.resize(cachetable.MMAX+1,cachetable.NORDER);

  forcetable.zforceC.resize(cachetable.MMAX+1,cachetable.NORDER);
  forcetable.zforceS.resize(cachetable.MMAX+1,cachetable.NORDER);

  if (z/cachetable.ASCALE > cachetable.Rtable) z =  cachetable.Rtable*cachetable.ASCALE;
  if (z/cachetable.ASCALE <-cachetable.Rtable) z = -cachetable.Rtable*cachetable.ASCALE;

  double X = (r_to_xi_cyl(r,cachetable.CMAPR,cachetable.ASCALE) - cachetable.XMIN)/cachetable.dX;
  double Y = (z_to_y_cyl(z,cachetable.CMAPZ,cachetable.HSCALE) - cachetable.YMIN)/cachetable.dY;

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

      forcetable.potC(mm,n) =
        (
        cachetable.potC[mm][n](ix  ,iy  ) * c00 +
         cachetable.potC[mm][n](ix+1,iy  ) * c10 +
         cachetable.potC[mm][n](ix  ,iy+1) * c01 +
         cachetable.potC[mm][n](ix+1,iy+1) * c11
        );

      forcetable.rforceC(mm,n) =
        (
        cachetable.rforceC[mm][n](ix  ,iy  ) * c00 +
        cachetable.rforceC[mm][n](ix+1,iy  ) * c10 +
        cachetable.rforceC[mm][n](ix  ,iy+1) * c01 +
        cachetable.rforceC[mm][n](ix+1,iy+1) * c11
        );

      forcetable.zforceC(mm,n) =
        (
        cachetable.zforceC[mm][n](ix  ,iy  ) * c00 +
        cachetable.zforceC[mm][n](ix+1,iy  ) * c10 +
        cachetable.zforceC[mm][n](ix  ,iy+1) * c01 +
        cachetable.zforceC[mm][n](ix+1,iy+1) * c11
        );

      // get sine values for m>0
      if (mm) {

  forcetable.potS(mm,n) =
    (
     cachetable.potS[mm][n](ix  ,iy  ) * c00 +
     cachetable.potS[mm][n](ix+1,iy  ) * c10 +
     cachetable.potS[mm][n](ix  ,iy+1) * c01 +
     cachetable.potS[mm][n](ix+1,iy+1) * c11
     );

        forcetable.rforceS(mm,n) =
    (
     cachetable.rforceS[mm][n](ix  ,iy  ) * c00 +
     cachetable.rforceS[mm][n](ix+1,iy  ) * c10 +
     cachetable.rforceS[mm][n](ix  ,iy+1) * c01 +
     cachetable.rforceS[mm][n](ix+1,iy+1) * c11
     );

        forcetable.zforceS(mm,n) =
    (
     cachetable.zforceS[mm][n](ix  ,iy  ) * c00 +
     cachetable.zforceS[mm][n](ix+1,iy  ) * c10 +
     cachetable.zforceS[mm][n](ix  ,iy+1) * c01 +
     cachetable.zforceS[mm][n](ix+1,iy+1) * c11
     );
      }

    } // end NORDER loop
  } // end MMAX loop

}




#endif
