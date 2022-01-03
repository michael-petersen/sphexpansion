/*
rotation matrices and helper rotation functions

MSP 30 Apr 2020 first version
MSP 14 May 2021 added header guard

what code does this mimic? we need to write down the convention somewhere...

 */
#ifndef ROTATE_H
#define ROTATE_H

#include <eigen3/Eigen/Dense>

using namespace Eigen;

//using namespace std;
using std::cout, std::cerr, std::endl;



MatrixXd return_euler_slater(double PHI, double THETA, double PSI, int BODY)
{
  double sph, cph, sth, cth, sps, cps;

  MatrixXd euler(3,3);

  sph = sin(PHI);
  cph = cos(PHI);
  sth = sin(THETA);
  cth = cos(THETA);
  sps = sin(PSI);
  cps = cos(PSI);
  
  euler(0,0) = -sps*sph + cth*cph*cps;
  euler(0,1) =  sps*cph + cth*sph*cps;
  euler(0,2) =  cps*sth;
      
  euler(1,0) = -cps*sph - cth*cph*sps;
  euler(1,1) =  cps*cph - cth*sph*sps;
  euler(1,2) = -sps*sth;
      
  euler(2,0) = -sth*cph;
  euler(2,1) = -sth*sph;
  euler(2,2) =  cth;

  if (BODY)
    return euler;//.Transpose();
  else
    return euler;
}

#endif
