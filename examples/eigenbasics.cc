#include <iostream>
#include <Eigen/Dense>

using Eigen::Tensor;

struct BasicMat
{
  Tensor m;
};

int main()
{
  BasicMat G;
  G.m.resize(2,2,2);
  G.m(0,0) = 3;
  G.m(1,0) = 2.5;
  G.m(0,1) = -1;
  //G.m(1,1) = m(1,0) + m(0,1);
  std::cout << G.m << std::endl;
}
