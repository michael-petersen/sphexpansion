/*
scaling.h

helper definitions for converting to scaled coordinates



 */
#ifndef SCALING_H
#define SCALING_H


#include <stdio.h>
#include <math.h>



double r_to_xi(double r, int cmap, double scale);

double xi_to_r(double xi, int cmap, double scale);

double d_xi_to_r(double xi, int cmap, double scale);

double z_to_y(double z, double hscale);

double y_to_z(double y, double hscale);




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


#endif
