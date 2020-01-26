#include "STABLE_NEO_HOOKEAN.h"
#include <iostream>

using namespace std;

STABLE_NEO_HOOKEAN::STABLE_NEO_HOOKEAN(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  _name = string("Stable-Neo-Hookean");
  std::cout << " Stable Neo-Hookean settings: " << std::endl;
  std::cout << " \t lambda: " << _lambda << std::endl;
  std::cout << " \t mu:     " << _mu << std::endl;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX STABLE_NEO_HOOKEAN::PK1(const MATRIX2& F)
{
  MATRIX2 DJ;
  DJ <<  F(1,1), -F(1,0),
        -F(0,1),  F(0,0);
  const Real J = F.determinant();
  const Real alpha = 1.0 + _mu / _lambda;

  MATRIX2 finalReturn = _mu * F + _lambda * (J - alpha) * DJ;
  return finalReturn;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX STABLE_NEO_HOOKEAN::PK2(const MATRIX2& F)
{
  return F.inverse() * PK1(F);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX STABLE_NEO_HOOKEAN::DPDF(const MATRIX& F)
{
  VECTOR f(4);
  f <<  F(1,1), -F(0,1),
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

  const Real alpha = 1.0 + _mu / _lambda;

  finalMatrix += _lambda * (f * f.transpose());
  finalMatrix += _lambda * (J - alpha) * anti;

  SelfAdjointEigenSolver<MATRIX> eigensolver(finalMatrix);
  VECTOR numerical = eigensolver.eigenvalues();
  MATRIX eigenvectors = eigensolver.eigenvectors();
  for (int x = 0; x < numerical.size(); x++)
    numerical[x] = (numerical[x] < 0.0) ? 0.0 : numerical[x];
  finalMatrix = eigenvectors * numerical.asDiagonal() * eigenvectors.transpose();

  return finalMatrix;
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR STABLE_NEO_HOOKEAN::flatten(const MATRIX& A)
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
Real STABLE_NEO_HOOKEAN::psi(const MATRIX2& F)
{
  MATRIX2 C = F.transpose() * F;
  Real Ic = C.trace();
  Real J = F.determinant();
  Real alpha = 1 + _mu / _lambda;

  return _mu * 0.5 * (Ic - 2) + _lambda * 0.5 * (J - alpha) * (J - alpha);
}
