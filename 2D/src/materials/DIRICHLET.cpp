#include "DIRICHLET.h"
#include <iostream>

using namespace std;

DIRICHLET::DIRICHLET(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  //_name = string("Neo-Hookean Length Preservation");
  //_name = string("L_Neo");
  _name = string("Dirichlet");
  std::cout << " Neo-Hookean length settings: " << std::endl;
  std::cout << " \t lambda: " << _lambda << std::endl;
  std::cout << " \t mu:     " << _mu << std::endl;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX DIRICHLET::PK1(const MATRIX2& F)
{
  return _mu * F;
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX DIRICHLET::PK2(const MATRIX2& F)
{
  MATRIX final(3,3);
  final.setZero();
  final(0,0) = _mu;
  final(1,1) = _mu;
  final(2,2) = _mu;
  return final;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX DIRICHLET::DPDF(const MATRIX& F)
{
  MATRIX final(4,4);
  int index = 0; 
  for (int j = 0; j < 2; j++)
    for (int i = 0; i < 2; i++, index++)
    {
      const MATRIX2 DFDFij = DFDF(F, i, j);
      MATRIX2 column = _mu * DFDFij;
      final.col(index) = flatten(column);
    }

  return final;
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR DIRICHLET::flatten(const MATRIX& A)
{
  VECTOR final(A.rows() * A.cols());

  int i = 0;
  for (int y = 0; y < A.cols(); y++)
    for (int x = 0; x < A.rows(); x++, i++)
      final[i] = A(x,y);

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
Real DIRICHLET::psi(const MATRIX2& F)
{
  MATRIX2 C = F.transpose() * F;
  Real Ic = C.trace();

  return _mu * 0.5 * (Ic - 3);
}
