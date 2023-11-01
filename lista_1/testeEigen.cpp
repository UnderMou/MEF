#include <iostream>
#include <Eigen/Dense>
 
using Eigen::MatrixXd;
using Eigen::VectorXd;
 
int main()
{
  MatrixXd m = MatrixXd::Random(2,2);
  m(0,0) = 10;
  m(0,1) = 0;
  m(1,0) = 0;
  m(1,1) = 10;

  //m = (m + MatrixXd::Constant(3,3,1.2)) * 50;
  std::cout << "m =" << std::endl << m << std::endl;
  VectorXd v(2);
  v << 1, 2;
  std::cout << "m * v =" << std::endl << m * v << std::endl;
}
