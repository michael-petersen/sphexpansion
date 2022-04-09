/*
clang++ --std=c++17 -I/opt/local/include -L/opt/local/lib -I../include/ -I/opt/local/include/eigen3 eigenbasics.cc -o obj/eigenbasics


*/
#include <iostream>
#include <Eigen/StdVector>
#include <Eigen/Dense>

using Eigen::MatrixXd;

struct BasicMat
{
  std::vector < std::vector< MatrixXd > > m;
};

int main()
{
  BasicMat G;
  G.m.resize(1);
  G.m[0].resize(1);
  G.m[0][0].resize(2,2);
  G.m[0][0](0,0) = 3;
  G.m[0][0](1,0) = 2.5;
  G.m[0][0](0,1) = -1;
  G.m[0][0](1,1) = G.m[0][0](1,0) + G.m[0][0](0,1);
  std::cout << G.m[0][0] << std::endl;


}
