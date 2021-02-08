/*
accumulate.h

compute coefficients from a distribution

MSP 30 Apr 2020 first commit

Plan:


 */

using namespace std;

// create 2- and 3-d array types
typedef boost::multi_array<double, 3> array_type3;
typedef boost::multi_array<double, 2> array_type2;

// create a spline array for the coefficients
typedef boost::multi_array<tk::spline, 2> spline_array;


struct AccCoefs
{
  int LMAX;            // the number of azimuthal harmonics
  int NMAX;            // the number of radial terms
  
  array_type2 coefs;   // the coefficient table, sized (LMAX+1)*(LMAX+1),NMAX

};


void accumulate_coefficients(array_type2 particles, SphCache cachetable, AccCoefs& coefmatrix) {
  /*
    accumulate coefficients

    pos is npart x 4, t, x, y, z

   */

  int l,m,loffset,moffset;
  
  double xx,yy,zz,mass,r2,r,costh,phi;
  array_type2 legs, dlegs;
  array_type2 potd,dpot;
  vector<double> cosm(cachetable.NMAX),sinm(cachetable.NMAX);
  
  int npart;

  double fac,fac1,fac2;//,fac4;

  double fac0=4.0*M_PI;
  
  int numl = (cachetable.LMAX+1)*(cachetable.LMAX+1);
  coefmatrix.coefs.resize(boost::extents[numl][cachetable.NMAX]);

  array_type2 normM;
  normM.resize(boost::extents[cachetable.LMAX+1][cachetable.NMAX]);
  for (int l=0; l<=cachetable.LMAX; l++) {	// with current binding from derived class
    for (int n=0; n<cachetable.NMAX; n++) {
      normM[l][n]  = 1.;//norm(n-1,l); // what is this exactly?? see SphericalBasis.cc
      // note: this is probably the spherical harmonic normalisation.
    }
  }
  
  
  npart = particles.shape()[0];
  
  // loop per particle
  for (int n=0;n<npart ; n++ ) {

    mass = particles[n][0];
    xx = particles[n][1];
    yy = particles[n][2];
    zz = particles[n][3];
    r2 = xx*xx + yy*yy + zz*zz;
    r = sqrt(r2+1.e-7); // this used to be a max call
    costh = zz/r;
    phi = atan2(yy,xx);

    dlegendre_R(cachetable.LMAX, costh, legs, dlegs);
    sinecosine_R(cachetable.LMAX, phi, cosm, sinm);
    get_dpotl(r, cachetable, potd, dpot);

    //fac = factorial[l][m] * legs[l][m]

    //fac4 = potd[l]*fac*fac0;

    // l loop
    for (l=0, loffset=0; l<=cachetable.LMAX; loffset+=(2*l+1), l++) {
	//		m loop
	for (m=0, moffset=0; m<=l; m++) {
	  if (m==0) {
	    for (n=0; n<cachetable.NMAX; n++) {
	      coefmatrix.coefs[loffset+moffset][n] += potd[l][n]*legs[l][m]*mass*fac0/normM[l][n];
	    }
	    moffset++;
	  }
	  else {
	    fac1 = legs[l][m]*cosm[m];
	    fac2 = legs[l][m]*sinm[m];

	    for (n=0; n<cachetable.NMAX; n++) {

	      coefmatrix.coefs[loffset+moffset  ][n] += potd[l][n]*fac1*mass*fac0/normM[l][n];
	      coefmatrix.coefs[loffset+moffset+1][n] += potd[l][n]*fac2*mass*fac0/normM[l][n];
	      
	    }
	    moffset+=2;
	  }
	} // end of m loop
    }// end of l loop

  } // end of particle loop
}



