/*
basis.h

Basic basis elements. See Basis.cc/.H for original implementation.

MSP 22 Apr 2020 clean version
MSP 14 May 2021 added header guard
MSP  9 Apr 2022 converted to Eigen

wishlist:
-update to featherstone & holmes algorithm

 */
#ifndef BASIS_H
#define BASIS_H

#include <Eigen/StdVector>
#include <Eigen/Dense>
using Eigen::MatrixXd;

// Machine constant for Legendre (note constexpr is not good in clang)
#ifndef EPSVAL
double EPS=std::numeric_limits<double>::min();
#define EPSVAL
#endif

int check_flags(int flag, int l)
{
  /*
  check binary flag for whether harmonic term is enabled

  monopole is always enabled!

  the rest of the values are a binary flag,
  2^0 = dipole
  2^1 = quadropole
  2^2 = octopole
  2^3 = l4
  etc

  for example, 111 -> 7 is 2^2 + 2^1 + 2^0, so l<=octopole are enabled.
               101 -> 5 is 2^2 +       2^0, so l=1,3 are enabled
               010 -> 2 is       2^1      , so l=2 only is enabled

  (in these examples, you set the integer value as flag, so 7, 5, 2 respectively.)

  */
  if (l==0) {
    return 1;
  } else {
    int lflag = int(pow(2,l-1));
    return (flag & lflag) == lflag;
  }
}


void legendre_R(int lmax, double x, MatrixXd& p)
{
  double fact, somx2, pll, pl1, pl2;
  int m, l;

  p(0,0) = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p(m,m) = pll;
      if (std::isnan(p(m,m))) std::cerr << "legendre_R: p[" << m << "][" << m << "]: pll=" << pll << "\n";
      fact += 2.0;
    }
  }

  for (m=0; m<lmax; m++) {
    pl2 = p(m,m);
    p(m+1,m) = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p(l,m) = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      if (std::isnan(p(l,m))) std::cerr << "legendre_R: p[" << l << "][" << m << "]: pll=" << pll << "\n";

      pl2 = pl1;
      pl1 = pll;
    }
  }

  if (std::isnan(x))
    std::cerr << "legendre_R: x\n";
  for(l=0; l<=lmax; l++)
    for (m=0; m<=l; m++)
      if (std::isnan(p(l,m))) std::cerr << "legendre_R: p[" << l << "][" << m << "] lmax=" << lmax << "\n";

}

void dlegendre_R(int lmax, double x, MatrixXd& p, MatrixXd& dp)
{
  double fact, somx2, pll, pl1, pl2;
  int m, l;

  p.resize(lmax+1,lmax+1);
  dp.resize(lmax+1,lmax+1);


  p(0,0) = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p(m,m) = pll;
      fact += 2.0;
    }
  }

  for (m=0; m<lmax; m++) {
    pl2 = p(m,m);
    p(m+1,m) = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p(l,m) = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      pl2 = pl1;
      pl1 = pll;
    }
  }

  if (1.0-fabs(x) < EPS) {
    if (x>0) x =   1.0 - EPS;
    else     x = -(1.0 - EPS);
  }

  somx2 = 1.0/(x*x - 1.0);
  dp(0,0) = 0.0;
  for (l=1; l<=lmax; l++) {
    for (m=0; m<l; m++)
      dp(l,m) = somx2*(x*l*p(l,m) - (l+m)*p(l-1,m));
    dp(l,l) = somx2*x*l*p(l,l);
  }
}

void sinecosine_R(int mmax, double phi, std::vector<double>& c, std::vector<double>& s)
{
  int m;

  c[0] = 1.0;
  s[0] = 0.0;

  c[1] = cos(phi);
  s[1] = sin(phi);

  for (m=2; m<=mmax; m++) {
    c[m] = 2.0*c[1]*c[m-1] - c[m-2];
    s[m] = 2.0*c[1]*s[m-1] - s[m-2];
  }
}

double factrl(int n)
{
    if(n > 1)
        return n * factrl(n - 1);
    else
        return 1;
}



void factorial (int lmax, Eigen::MatrixXd& factorial) {

  if (lmax>1) {

    factorial.resize(lmax+1,lmax+1);

    for (int l=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++)
        factorial(l,m) = factrl(l-m)/factrl(l+m);
    }
  } else {

  }

}

#endif
