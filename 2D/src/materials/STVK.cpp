#include "STVK.h"
#include <iostream>

STVK::STVK(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  //_lambda = 0;

  _name = std::string("StVK");
  /*
  std::cout << " StVK settings: " << std::endl;
  std::cout << " \t lambda: " << _lambda << std::endl;
  std::cout << " \t mu:     " << _mu << std::endl;
  */
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX STVK::PK1(const MATRIX2& F)
{
  return F * PK2(F);
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX STVK::PK2(const MATRIX2& F)
{
  MATRIX2 E = 0.5 * (F.transpose() * F - MATRIX2::Identity());
  MATRIX2 S = _lambda * E.trace() * MATRIX2::Identity() + 2.0 * _mu * E;
  return S;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of S, the second Piola-Kirchhoff
// with respect to E = 1/2 (F^T * F - I)
///////////////////////////////////////////////////////////////////////
MATRIX STVK::DSDE()
{
  MATRIX final(4,4);
  final.setZero();

  final(0,0) = _lambda + 2.0 * _mu;
  final(3,0) = _lambda;

  final(1,1) = 2.0 * _mu;

  final(2,2) = 2.0 * _mu;

  final(0,3) = _lambda;
  final(3,3) = _lambda + 2.0 * _mu;

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the gradient of FS (i.e. PK1) w.r.t F assuming that S is frozen
///////////////////////////////////////////////////////////////////////
MATRIX STVK::DFSDF(const MATRIX& S)
{
  MATRIX final(4,4);
  final.setZero();

  final(0,0) = S(0,0);
  final(1,1) = S(0,0);
  final(2,2) = S(1,1);
  final(3,3) = S(1,1);

  final(2,0) = S(0,1);
  final(3,1) = S(0,1);
  
  final(0,2) = S(1,0);
  final(1,3) = S(1,0);

  return final;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX STVK::DPDF(const MATRIX& F)
{
  MATRIX S = PK2(F);

  MATRIX IS = DFSDF(S);

  MATRIX diagF = blockDiag(F, 2);

  // no getting around computing this one -- changes with F
  MATRIX DEDF = computeDEDF(F);

  return IS + diagF * DSDE() * DEDF;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of E = 1/2 (F^T * F - I)
// with respect to F
///////////////////////////////////////////////////////////////////////
MATRIX STVK::computeDEDF(const MATRIX2& F)
{
  const Real& f00 = F(0,0);
  const Real& f10 = F(1,0);
  const Real& f01 = F(0,1);
  const Real& f11 = F(1,1);
  MATRIX final(4,4);

  final(0,0) = 2.0 * f00;
  final(1,0) = f01;
  final(2,0) = f01;
  final(3,0) = 0.0;

  final(0,1) = 2.0 * f10;
  final(1,1) = f11;
  final(2,1) = f11;
  final(3,1) = 0.0;

  final(0,2) = 0.0;
  final(1,2) = f00;
  final(2,2) = f00;
  final(3,2) = 2.0 * f01;
  
  final(0,3) = 0.0;
  final(1,3) = f10;
  final(2,3) = f10;
  final(3,3) = 2.0 * f11;

  // E = 0.5 (F^T * F - I), so a 0.5 factor is needed in front
  final *= 0.5;

  return final;
}

///////////////////////////////////////////////////////////////////////
// repeat the given matrix "repeat" times along the block diagonal
///////////////////////////////////////////////////////////////////////
MATRIX STVK::blockDiag(const MATRIX& A, const int repeats)
{
  const int rows = A.rows();
  const int cols = A.cols();
  MATRIX final(rows * repeats, cols * repeats);
  final.setZero();

  for (int y = 0; y < cols; y++)
    for (int x = 0; x < rows; x++)
    {
      for (int i = 0; i < repeats; i++)
        final(i * rows + x, i * cols + y) = A(x,y);
    }

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
Real STVK::psi(const MATRIX2& F)
{
  MATRIX2 E = 0.5 * (F.transpose() * F - MATRIX2::Identity());
  return _mu * E.squaredNorm() + _lambda * 0.5 * pow(E.trace(), 2.0); 
  //return _mu * E.squaredNorm() + _lambda * 0.5 * E.trace();
  //return _lambda * 0.5 * pow(E.trace(), 2.0);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real STVK::constant() const
{
  //return -_lambda * 0.5;
  //return 4;
  return _lambda / 2;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real STVK::length(const Real sigma) const 
{
  return _mu * 0.25 * pow(sigma * sigma - 1.0, 2.0) +
         _lambda / 8.0 * (pow(sigma, 4.0) - 4.0 * sigma * sigma);
         //_lambda * 0.25 * sigma * sigma;
         //_lambda * 0.25 * (sigma * sigma); 
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real STVK::area(const Real s1s2) const
{
  return _lambda * 0.25 * s1s2 * s1s2;
  //return 0;
}
