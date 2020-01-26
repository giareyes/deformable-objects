#include "LENGTH_STVK.h"
#include <iostream>

using namespace std;

LENGTH_STVK::LENGTH_STVK(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  //_name = string("STVK Length Preservation");
  _name = string("L_StVK");
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = 2 mu F E
///////////////////////////////////////////////////////////////////////
MATRIX LENGTH_STVK::PK1(const MATRIX2& F)
{
  MATRIX2 E = 0.5 * (F.transpose() * F - MATRIX2::Identity());
  return 2.0 * _mu * F * E;
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor = 2 mu E
///////////////////////////////////////////////////////////////////////
MATRIX LENGTH_STVK::PK2(const MATRIX2& F)
{
  MATRIX2 E = 0.5 * (F.transpose() * F - MATRIX2::Identity());
  return 2.0 * _mu * E;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX LENGTH_STVK::DPDF(const MATRIX& F)
{
  MATRIX2 E = 0.5 * (F.transpose() * F - MATRIX2::Identity());

  MATRIX final(4,4);
  int index = 0; 
  for (int j = 0; j < 2; j++)
    for (int i = 0; i < 2; i++, index++)
    {
      const MATRIX2 DFDFij = DFDF(F, i, j);
      const MATRIX2 DEDFij = DEDF(F, i, j);
      MATRIX2 column = 2.0 * _mu * (DFDFij * E + F * DEDFij);
      final.col(index) = flatten(column);
    }

  return final;
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR LENGTH_STVK::flatten(const MATRIX& A)
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
Real LENGTH_STVK::psi(const MATRIX2& F)
{
  MATRIX2 E = 0.5 * (F.transpose() * F - MATRIX2::Identity());
  return _mu * E.squaredNorm();
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of E = 1/2 (F^T * F - I)
// with respect to F
///////////////////////////////////////////////////////////////////////
MATRIX2 LENGTH_STVK::DEDF(const MATRIX2& F, int i, int j)
{
  const Real& f00 = F(0,0);
  const Real& f10 = F(1,0);
  const Real& f01 = F(0,1);
  const Real& f11 = F(1,1);
  MATRIX2 final;

  if (i == 0 && j == 0)
  {
    final(0,0) = 2.0 * f00;
    final(1,0) = f01;
    final(0,1) = f01;
    final(1,1) = 0.0;
  }
  if (i == 1 && j == 0)
  {
    final(0,0) = 2.0 * f10;
    final(1,0) = f11;
    final(0,1) = f11;
    final(1,1) = 0.0;
  }

  if (i == 0 && j == 1)
  {
    final(0,0) = 0.0;
    final(1,0) = f00;
    final(0,1) = f00;
    final(1,1) = 2.0 * f01;
  }
  
  if (i == 1 && j == 1)
  {
    final(0,0) = 0.0;
    final(1,0) = f10;
    final(0,1) = f10;
    final(1,1) = 2.0 * f11;
  }
  return 0.5 * final;
}
