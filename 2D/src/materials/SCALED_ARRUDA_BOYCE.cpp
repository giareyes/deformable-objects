#include "SCALED_ARRUDA_BOYCE.h"
#include <iostream>

using namespace std;

SCALED_ARRUDA_BOYCE::SCALED_ARRUDA_BOYCE(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  _lambda = 1.0;
  //_lambda = 12.1;
  _mu = 2.0;
  _beta = 1.0;

  //_alpha = 0.5 + sqrt(_lambda * (_lambda - 4.0 * _mu - 8.0 * _beta)) / (2.0 * _lambda);
  //_alpha = sqrt((_lambda - _mu) / (2.0 * _beta + _lambda));
  _alpha = 1.0;
  Real firstAlpha = 0.5 + sqrt(_lambda * (_lambda - 4.0 * _mu - 8.0 * _beta)) / (2.0 * _lambda);
  Real secondAlpha = 0.5 - sqrt(_lambda * (_lambda - 4.0 * _mu - 8.0 * _beta)) / (2.0 * _lambda);

  _name = string("Scaled-Arruda-Boyce");
  std::cout << " Scaled Arruda-Boyce settings: " << std::endl;
  std::cout << " \t lambda: " << _lambda << std::endl;
  std::cout << " \t mu:     " << _mu << std::endl;
  std::cout << " \t alpha:  " << _alpha << std::endl;
  std::cout << " \t beta:   " << _beta << std::endl;
  std::cout << " \t first alpha:  " << firstAlpha << std::endl;
  std::cout << " \t second alpha: " << secondAlpha << std::endl;

  Real stability = _mu / (_lambda * _alpha) - 1.0;
  std::cout << "\t stability param: " << stability << endl;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX SCALED_ARRUDA_BOYCE::PK1(const MATRIX2& F)
{
  /*
  MATRIX2 DJ;
  DJ <<  F(1,1), -F(1,0),
        -F(0,1),  F(0,0);
  const Real J = F.determinant();
  const MATRIX2 C = F.transpose() * F;
  const Real Ic = C.trace();

  MATRIX2 finalReturn = _mu * F + + _beta * Ic * F + _lambda * _alpha * (_alpha * J - 1.0) * DJ;
  return finalReturn;
  */
  MATRIX2 DJ;
  DJ <<  F(1,1), -F(1,0),
        -F(0,1),  F(0,0);
  const Real J = F.determinant();
  const MATRIX2 C = F.transpose() * F;
  const Real Ic = C.trace();

  MATRIX2 finalReturn = _mu * _alpha * _alpha * F + 
                        _beta * pow(_alpha, 4.0) * Ic * F + 
                        _lambda * _alpha * _alpha * (_alpha * _alpha * J - 1.0) * DJ;
  return finalReturn;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX SCALED_ARRUDA_BOYCE::PK2(const MATRIX2& F)
{
  return F.inverse() * PK1(F);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX SCALED_ARRUDA_BOYCE::DPDF(const MATRIX& F)
{
  VECTOR g(4);
  g <<  F(1,1), -F(0,1),
        -F(1,0),  F(0,0);
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
  Real alphaSq = _alpha * _alpha;
  Real alphaQuar = alphaSq * alphaSq;
  for (int x = 0; x < 4; x++)
    finalMatrix(x,x) = _mu * alphaSq + _beta * Ic * alphaQuar;

  finalMatrix += _lambda * alphaQuar * (g * g.transpose());
  finalMatrix += _lambda * alphaSq * (alphaSq * J - 1.0) * anti;
  
  VECTOR f = flatten(F);
  finalMatrix += 2.0 * _beta * alphaQuar * (f * f.transpose());

  // peek at the eigenvalues
  VECTOR eigenvalues = finalMatrix.eigenvalues().real();
  bool hasNegative = false;
  for (int x = 0; x < 4; x++)
    if (eigenvalues[x] < -1e-8)
      hasNegative = true;

  if (hasNegative)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Negative found: " << eigenvalues.transpose() << endl;
    exit(0);
  }

  return finalMatrix;
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR SCALED_ARRUDA_BOYCE::flatten(const MATRIX& A)
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
Real SCALED_ARRUDA_BOYCE::psi(const MATRIX2& F)
{
  MATRIX2 C = F.transpose() * F;
  Real Ic = C.trace();
  Real J = F.determinant();

  /*
  return _mu * 0.5 * (Ic - 3) + 
         _beta * 0.25 * (Ic * Ic - 9) +
         _lambda * 0.5 * (_alpha * J - 1.0) * (_alpha * J - 1.0);
         */

  return _mu * 0.5 * (_alpha * _alpha) * (Ic - 3) + 
         _beta * 0.25 * (pow(_alpha, 4.0) * Ic * Ic - 9) +
         _lambda * 0.5 * (_alpha * _alpha * J - 1.0) * (_alpha * _alpha * J - 1.0);
}
