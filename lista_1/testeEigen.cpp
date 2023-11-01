#include <iostream>
#include <Eigen/Dense>
 
using Eigen::MatrixXd;
using Eigen::VectorXd;
 
int main()
{
  MatrixXd A = MatrixXd::Random(3,3);
  A(0,0) = 1;
  A(0,1) = 2;
  A(0,2) = 3;
  A(1,0) = 4;
  A(1,1) = 5;
  A(1,2) = 6;
  A(2,0) = 7;
  A(2,1) = 8;
  A(2,2) = 10;

  VectorXd v(3);
  v << 3, 3, 4;
  
  std::cout << "Here is the matrix A:\n" << A << std::endl;
  std::cout << "Here is the vector b:\n" << v << std::endl;
  Eigen::Vector3f x = A.colPivHouseholderQr().solve(v);
  std::cout << "The solution is:\n" << x << std::endl;
}

