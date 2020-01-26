#include "SCALED_NEO_HOOKEAN.h"
#include <iostream>

using namespace std;

SCALED_NEO_HOOKEAN::SCALED_NEO_HOOKEAN(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  _alpha = 0.5 + sqrt(_lambda - 4.0 * _mu) / (2.0 * sqrt(_lambda));
  Real otherAlpha = 0.5 - sqrt(_lambda - 4.0 * _mu) / (2.0 * sqrt(_lambda));

  _name = string("Scaled-Neo-Hookean");
  std::cout << " Scaled Neo-Hookean settings: " << std::endl;
  std::cout << " \t lambda: " << _lambda << std::endl;
  std::cout << " \t mu:     " << _mu << std::endl;
  std::cout << " \t alpha:  " << _alpha << std::endl;
  std::cout << " \t other alpha:  " << otherAlpha << std::endl;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX SCALED_NEO_HOOKEAN::PK1(const MATRIX2& F)
{
  MATRIX2 DJ;
  DJ <<  F(1,1), -F(1,0),
        -F(0,1),  F(0,0);
  const Real J = F.determinant();

  MATRIX2 finalReturn = _mu * F + _lambda * _alpha * (_alpha * J - 1.0) * DJ;
  return finalReturn;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX SCALED_NEO_HOOKEAN::PK2(const MATRIX2& F)
{
  return F.inverse() * PK1(F);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX SCALED_NEO_HOOKEAN::DPDF(const MATRIX& F)
{
  VECTOR g(4);
  g <<  F(1,1), -F(0,1),
        -F(1,0),  F(0,0);
  const Real J = F.determinant();

  MATRIX anti(4,4);
  anti.setZero();
  anti(0,3) = 1;
  anti(1,2) = -1;
  anti(2,1) = -1;
  anti(3,0) = 1;

  MATRIX finalMatrix(4,4);
  finalMatrix.setZero();
  for (int x = 0; x < 4; x++)
    finalMatrix(x,x) = _mu;

  finalMatrix += _lambda * _alpha * _alpha * (g * g.transpose());
  finalMatrix += _lambda * _alpha * (_alpha * J - 1.0) * anti;

  return finalMatrix;
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR SCALED_NEO_HOOKEAN::flatten(const MATRIX& A)
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
Real SCALED_NEO_HOOKEAN::psi(const MATRIX2& F)
{
  MATRIX2 C = F.transpose() * F;
  Real Ic = C.trace();
  Real J = F.determinant();

  return _mu * 0.5 * (Ic - 3) + _lambda * 0.5 * (_alpha * J - 1.0) * (_alpha * J - 1.0);
}
