/*
scaling.h

helper definitions for converting to scaled coordinates

-check these definitions
-convert r_to_xi calls to r_to_xi_sph calls

 */
#ifndef SCALING_H
#define SCALING_H


#include <stdio.h>
#include <math.h>


// spherical scalings
double r_to_xi(double r, int cmap, double scale);
double xi_to_r(double xi, int cmap, double scale);
double d_xi_to_r(double xi, int cmap, double scale);

// cylindrical scalings
double r_to_xi_cyl  (double r,  int cmapr, double ascale);
double xi_to_r_cyl  (double xi, int cmapr, double ascale);
double d_xi_to_r_cyl(double xi, int cmapr, double ascale);
double z_to_y_cyl   (double y,  int cmapz, double hscale);
double y_to_z_cyl   (double y,  int cmapz, double hscale);
double d_y_to_z_cyl (double y,  int cmapz, double hscale);
double z_to_y(double z, double hscale);
double y_to_z(double y, double hscale);


double r_to_xi_cyl(double r, int cmapr, double ascale)
{
  if (cmapr>0) {
    if (r<0.0) {
      std::cerr << "Issue in r mapping..." << std::endl;
    }
    return (r/ascale - 1.0)/(r/ascale + 1.0);
  } else {
    if (r<0.0)  {
      std::cerr << "Issue in r mapping..." << std::endl;
    }
    return r;
  }
}

double xi_to_r_cyl(double xi, int cmapr, double ascale)
{
  if (cmapr>0) {
    if (xi<-1.0) std::cerr << "Issue in xi mapping..." << std::endl;
    if (xi>=1.0) std::cerr << "Issue in xi mapping..." << std::endl;

    return (1.0 + xi)/(1.0 - xi) * ascale;
  } else {
    return xi;
  }

}


double d_xi_to_r_cyl(double xi, int cmapr, double ascale)
{
  if (cmapr>0) {
    if (xi<-1.0) std::cerr << "Issue in xi mapping..." << std::endl;
    if (xi>=1.0) std::cerr << "Issue in xi mapping..." << std::endl;

    return 0.5*(1.0 - xi)*(1.0 - xi) / ascale;
  } else {
    return 1.0;
  }
}


double r_to_xi(double r, int cmap, double scale)
{
  if ( cmap == 1 ) {

    if (r<0.0) {
      printf("radius < 0!");
      return 0.0;
	} else {
      return (r/scale-1.0)/(r/scale+1.0);
    }

  } else if ( cmap == 2 ) {
    if (r<=0.0) {
      printf("radius <= 0!");
      return 0.0;
    } else {
    return log(r);
    }

  } else {

    if (r<0.0) {
      printf("radius < 0!");
      return 0.0;
    } else {

    return r;
    }
  }
}



double xi_to_r(double xi, int cmap, double scale)
{
  if ( cmap == 1 ) {
    if (xi<-1.0) printf("xi < -1!");
    if (xi>=1.0) printf("xi >= 1!");

    return (1.0+xi)/(1.0 - xi) * scale;

     } else if (cmap==2) {
    return exp(xi);

  } else {
    return xi;
  }

}


double d_xi_to_r(double xi, int cmap, double scale)
{
  if ( cmap == 1 ) {
    if (xi<-1.0) printf("xi < -1!");
    if (xi>=1.0) printf("xi >= 1!");

    return 0.5*(1.0-xi)*(1.0-xi)/scale;

      } else if (cmap==2) {
    return exp(-xi);

  } else {
    return 1.0;
  }
}


double z_to_y(double z, double hscale)
{
  return z /( fabs(z)+1.e-10) * asinh( fabs(z/hscale));
}

double y_to_z(double y, double hscale)
{
  return hscale*sinh(y);
}

double z_to_y_cyl(double z, int cmapz, double hscale)
{
  if (cmapz==1)
    return z/(fabs(z)+std::numeric_limits<double>::min())*asinh(fabs(z/hscale));
  else if (cmapz==2)
    return z/sqrt(z*z + hscale*hscale);
  else
    return z;
}

double y_to_z_cyl(double y, int cmapz, double hscale)
{
  // Compute Z from non-dimensional vertical coordinate
  if (cmapz==1)
    return hscale*sinh(y);
  else if (cmapz==2) {
    return y * hscale/sqrt(1.0 - y*y);
  } else
    return y;
}


double d_y_to_z_cyl(double y, int cmapz, double hscale)
{
  // For measure transformation in vertical coordinate
  if (cmapz==1)
    return hscale*cosh(y);
  else if (cmapz==2) {
    return hscale*pow(1.0-y*y, -1.5);
  } else
    return 1.0;
}

#endif
