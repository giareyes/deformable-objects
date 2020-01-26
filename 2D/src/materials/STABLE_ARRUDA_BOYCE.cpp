#include "STABLE_ARRUDA_BOYCE.h"
#include <iostream>

using namespace std;

STABLE_ARRUDA_BOYCE::STABLE_ARRUDA_BOYCE(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu), _beta(1.0)
{
  _name = string("Stable-Arruda-Boyce");
  std::cout << " Stable Arruda-Boyce settings: " << std::endl;

  //_lambda = 1.0;
  _lambda = 1.0;
  _mu = 1.0;
  _beta = 0.25;

  //_alpha = 1.0 + (_mu / _lambda) * (1.0 + 2.0 * _beta);
  _alpha = 1.0;

  std::cout << " \t lambda: " << _lambda << std::endl;
  std::cout << " \t mu:     " << _mu << std::endl;
  std::cout << " \t beta:   " << _beta << std::endl;
  std::cout << " \t alpha:  " << _alpha<< std::endl;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX STABLE_ARRUDA_BOYCE::PK1(const MATRIX2& F)
{
  MATRIX2 DJ;
  DJ <<  F(1,1), -F(1,0),
        -F(0,1),  F(0,0);
  const Real J = F.determinant();
  const MATRIX2 C = F.transpose() * F;
  const Real Ic = C.trace();

  MATRIX2 finalReturn = (_mu + _beta * Ic) * F + _lambda * (J - _alpha) * DJ;
  return finalReturn;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX STABLE_ARRUDA_BOYCE::PK2(const MATRIX2& F)
{
  return F.inverse() * PK1(F);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX STABLE_ARRUDA_BOYCE::DPDF(const MATRIX& F)
{
  VECTOR g(4);
  g <<  F(1,1), -F(0,1),
        -F(1,0),  F(0,0);

  VECTOR f = flatten(F);
  const Real J = F.determinant();
  const MATRIX2 C = F.transpose() * F;
  const Real Ic = C.trace();

  MATRIX anti(4,4);
  anti.setZero();
  anti(0,3) = 1;
  anti(1,2) = -1;
  anti(2,1) = -1;
  anti(3,0) = 1;

  MATRIX finalMatrix(4,4);
  finalMatrix.setZero();
  for (int x = 0; x < 4; x++)
    finalMatrix(x,x) = (_mu + _beta * Ic);

  finalMatrix += _lambda * (g * g.transpose());
  finalMatrix += _lambda * (J - _alpha) * anti;
  finalMatrix += 2.0 * _beta * (f * f.transpose());

  // peek at the eigenvalues
  VECTOR eigenvalues = finalMatrix.eigenvalues().real();
  bool hasNegative = false;
  for (int x = 0; x < 4; x++)
    if (eigenvalues[x] < 0.0)
      hasNegative = true;

  if (hasNegative)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Negative found: " << eigenvalues.transpose() << endl;
    //exit(0);
  }

  return finalMatrix;
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR STABLE_ARRUDA_BOYCE::flatten(const MATRIX& A)
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
Real STABLE_ARRUDA_BOYCE::psi(const MATRIX2& F)
{
  MATRIX2 C = F.transpose() * F;
  Real Ic = C.trace();
  Real J = F.determinant();

  return _mu * 0.5 * (Ic - 3) + (_beta * 0.25) * (Ic * Ic - 9) + _lambda * 0.5 * (J - _alpha) * (J - _alpha);
}
