/*
transform.h

helper functions driving transformations from (x,y,z) to (rho,phi,theta)

MSP 22 Apr 2020 clean version
MSP 24 Apr 2020 revised to handle edge cases 
MSP 14 May 2021 added header guard


wishlist:
-vector versions of the functions to stop having to pass so much text

 */
#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <math.h>

double EPS=1.e-10;


void cartesian_to_cylindrical(double  x, double    y,
			      double& r, double& phi)
{
  r     = sqrt(x*x + y*y);
  phi   = atan2(y,x);
}


void cylindrical_to_cartesian(double  r, double  phi,
			      double& x, double&   y)
{
  x = r * cos(phi);
  y = r * sin(phi);
}

void cylindrical_forces_to_cartesian(double  r, double phi,
				     double fr, double fp,
	                             double& fx, double& fy)
{

  double x,y;
  cylindrical_to_cartesian(r, phi, x, y);

  // check the handedness of the coordinate system.
  
  fx = (x*fr - y*fp)/r;

  fy = (y*fr + x*fp)/r;

}



void cartesian_to_spherical(double  x, double    y, double      z,
			    double& r, double& phi, double& theta)
{
  r     = sqrt(x*x + y*y + z*z);
  phi   = atan2(y,x);
  theta = acos( z/r );  
}

void spherical_to_cartesian(double  r, double  phi, double  theta,
			    double& x, double&   y, double&     z)
{
  x = r * sin(theta) * cos(phi);
  y = r * sin(theta) * sin(phi);
  z = r * cos(theta);
}

void spherical_forces_to_cartesian_legacy(double r3, double phi, double theta,
				   double fr, double fp, double ft,
	                           double& fx, double& fy, double& fz)
{

  double x,y,z;
  spherical_to_cartesian(r3, phi, theta, x, y, z);
  double r2 = sqrt(x*x + y*y);

  fx = ( x*(r2*fr + z*ft) - y*r3*fp )/(r2*r3);

  fy = ( y*(r2*fr + z*ft) + x*r3*fp )/(r2*r3);
    
  fz = ( z*fr - r2*ft )/(r3);
  
}


void spherical_forces_to_cartesian(double r3, double phi, double theta,
				   double fr, double fp, double ft,
	                           double& fx, double& fy, double& fz)
{

  double x,y,z;
  spherical_to_cartesian(r3, phi, theta, x, y, z);
  double r2 = sqrt(x*x + y*y);

  // check the handedness of the coordinate system.
  
  fx = - (( fr*(x/r3) - ft*(x*z/(r3*r3*r3))) + fp*(y/(r2*r2)));

  fy = - (( fr*(y/r3) - ft*(y*z/(r3*r3*r3))) - fp*(x/(r2*r2)));
    
  fz = - ( fr*(z/r3) + ft*( (r2*r2)/(r3*r3*r3)) );
  
}

#endif
