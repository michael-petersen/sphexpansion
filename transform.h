/*
transform.h

helper functions driving transformations from (x,y,z) to (rho,phi,theta)

clean version, MSP 22 April 2020


 */

#include <math.h>

void cartesian_to_spherical(double  x, double    y, double      z,
			    double& r, double& phi, double& theta)
{
  r     = sqrt(x*x + y*y + z*z);
  phi   = atan2(y,x);
  theta = atan(sqrt(x*x + y*y)/z);  
}

void spherical_to_cartesian(double  r, double  phi, double  theta,
			    double& x, double&   y, double&     z)
{
  if (theta<1.e-6) theta = 1.e-6;
  x = r * sin(theta) * cos(phi);
  y = r * sin(theta) * sin(phi);
  z = r * cos(theta);
}

void spherical_forces_to_cartesian(double r3, double phi, double theta,
				   double fr, double fp, double ft,
	                           double& fx, double& fy, double& fz)
{

  double x,y,z;
  spherical_to_cartesian(r3, phi, theta, x, y, z);
  double r2 = sqrt(x*x + y*y);

  
  fx = ( fr*(x/r3) - ft*(x*z/(r3*r3*r3))) + fp*(y/(r2*r2));

  fy = ( fr*(y/r3) - ft*(y*z/(r3*r3*r3))) - fp*(x/(r2*r2));
    
  fz = ( fr*(z/r3) + ft*( (r2*r2)/(r3*r3*r3)) );
  
}

