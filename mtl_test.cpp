/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

#include "IdentityMatrix.hpp"
//#include "GraphSymmetricMatrix.hpp"
int main()
{
  const int size=50;
  IdentityMatrix A= IdentityMatrix(size);
  mtl::dense_vector<double> x(size,0.0);
  mtl::dense_vector<double> b(size,1.0);
  assert(x!=b);
  itl::pc::identity<IdentityMatrix> P(A);
  itl::noisy_iteration<double > iter(b, 500, 1.e-6);
  itl::cg(A,x,b,P,iter);
  assert(x==b);
  return 0;
}
