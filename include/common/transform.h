/*
transform.h

helper functions driving transformations from (x,y,z) to (rho,phi,theta)
  all functions have overloaded ArrayXd versions; can clean up elsewhere in code base

MSP 22 Apr 2020 clean version
MSP 24 Apr 2020 revised to handle edge cases
MSP 14 May 2021 added header guard
MSP 18 Apr 2022 add overloaded versions using Eigen

wishlist:
-add coordinate handedness options

 */
#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <math.h>
#include <cmath>

// Eigen VectorXd
#include <Eigen/Dense>
using Eigen::ArrayXd;

#ifndef EPSVAL
double EPS=std::numeric_limits<double>::min();
#define EPSVAL
#endif



void cartesian_to_cylindrical(double  x, double    y,
			      									double& r, double& phi)
{
  r     = sqrt(x*x + y*y);
  phi   = atan2(y,x);
}

void cartesian_to_cylindrical(ArrayXd  x, ArrayXd    y,
															ArrayXd& r, ArrayXd& phi)
{
	// OVERLOADED
  r     = sqrt(x*x + y*y);

	// no atan2 in eigen, so we have to loop.
	for (int i=0; i<phi.size(); i++) {
  	phi(i)   = std::atan2(y(i),x(i));
	}
}


void cylindrical_to_cartesian(double  r, double  phi,
												      double& x, double&   y)
{
  x = r * cos(phi);
  y = r * sin(phi);
}

void cylindrical_to_cartesian(ArrayXd  r, ArrayXd  phi,
												      ArrayXd& x, ArrayXd&   y)
{
	// OVERLOADED
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

void cylindrical_forces_to_cartesian(ArrayXd  r, ArrayXd phi,
				     												 ArrayXd fr, ArrayXd fp,
	                             			 ArrayXd& fx,ArrayXd& fy)
{
  // OVERLOADED
  ArrayXd x,y;
  cylindrical_to_cartesian(r, phi, x, y);

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

void cartesian_to_spherical(ArrayXd  x, ArrayXd    y, ArrayXd      z,
			    									ArrayXd& r, ArrayXd& phi, ArrayXd& theta)
{
	// OVERLOADED
  r     = sqrt(x*x + y*y + z*z);

	// no atan2 in eigen, so we have to loop.
	for (int i=0; i<phi.size(); i++) {
    phi(i)   = atan2(y(i),x(i));
    theta(i) = acos( z(i)/r(i) );
	}
}

void spherical_to_cartesian(double  r, double  phi, double  theta,
			    									double& x, double&   y, double&     z)
{
  x = r * sin(theta) * cos(phi);
  y = r * sin(theta) * sin(phi);
  z = r * cos(theta);
}


void spherical_to_cartesian(ArrayXd  r, ArrayXd  phi, ArrayXd  theta,
			    									ArrayXd& x, ArrayXd&   y, ArrayXd&     z)
{
	// OVERLOADED
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

void spherical_forces_to_cartesian_legacy(ArrayXd r3, ArrayXd phi, ArrayXd theta,
				   																ArrayXd fr, ArrayXd fp,  ArrayXd ft,
	                           							ArrayXd& fx,ArrayXd& fy, ArrayXd& fz)
{
	// OVERLOADED

  ArrayXd x,y,z;
  spherical_to_cartesian(r3, phi, theta, x, y, z);
  ArrayXd r2 = sqrt(x*x + y*y);

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

  fx = - (( fr*(x/r3) - ft*(x*z/(r3*r3*r3))) + fp*(y/(r2*r2)));
  fy = - (( fr*(y/r3) - ft*(y*z/(r3*r3*r3))) - fp*(x/(r2*r2)));
  fz = - (  fr*(z/r3) + ft*( (r2*r2)/(r3*r3*r3)) );

}

void spherical_forces_to_cartesian(ArrayXd r3, ArrayXd phi, ArrayXd theta,
				   												 ArrayXd fr, ArrayXd fp,  ArrayXd ft,
	                           			 ArrayXd& fx,ArrayXd& fy, ArrayXd& fz)
{
	// OVERLOADED
  ArrayXd x,y,z;
  spherical_to_cartesian(r3, phi, theta, x, y, z);
  ArrayXd r2 = sqrt(x*x + y*y);

  fx = - (( fr*(x/r3) - ft*(x*z/(r3*r3*r3))) + fp*(y/(r2*r2)));
  fy = - (( fr*(y/r3) - ft*(y*z/(r3*r3*r3))) - fp*(x/(r2*r2)));
  fz = - (  fr*(z/r3) + ft*( (r2*r2)/(r3*r3*r3)) );

}

#endif
