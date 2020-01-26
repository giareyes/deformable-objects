#include "VOLUME_MR.h"
#include <iostream>

using namespace std;

VOLUME_MR::VOLUME_MR(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  _ratio = _mu / _lambda;

  //_name = string("Mooney-Rivlin Volume Preservation");
  _name = string("V_Moon");
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX VOLUME_MR::PK1(const MATRIX2& F)
{
  const Real J = F.determinant();
  MATRIX2 DJDF;
  DJDF <<  F(1,1), -F(1,0),
          -F(0,1), F(0,0);

  //return _lambda * (J - 1) * DJDF; 
  return _lambda * (J - 1 - _ratio) * DJDF; 
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX VOLUME_MR::PK2(const MATRIX2& F)
{
  return F.inverse() * PK1(F);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX VOLUME_MR::DPDF(const MATRIX& F)
{
  MATRIX final(4,4);
  int index = 0; 
  for (int j = 0; j < 2; j++)
    for (int i = 0; i < 2; i++, index++)
    {
      MATRIX column = DPvDF(F, i, j);
      final.col(index) = flatten(column);
    }
  return final;
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR VOLUME_MR::flatten(const MATRIX& A)
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
Real VOLUME_MR::psi(const MATRIX2& F)
{
  Real J = F.determinant();
  //return _lambda * 0.5 * (J - 1) * (J - 1);
  return _lambda * 0.5 * (J - 1 - _ratio) * (J - 1 - _ratio);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX2 VOLUME_MR::DPvDF(const MATRIX2& F, int i, int j)
{
  MATRIX2 A;
  A <<  F(1,1), -F(1,0),
       -F(0,1),  F(0,0);
  
  Real J = F.determinant();
  MATRIX2 final;

  MATRIX2 zeros;
  zeros.setZero();

  Real coeff; 
  if (i == 0)
  {
    // F(0,0)
    if (j == 0)
    {
      coeff = F(1,1);
      zeros(1,1) = 1;
    }
    // F(0,1)
    else
    {
      coeff = -F(1,0);
      zeros(1,0) = -1;
    }
  }
  else
  {
    // F(1,0)
    if (j == 0)
    {
      coeff = -F(0,1);
      zeros(0,1) = -1;
    }
    // F(1,1)
    else
    {
      coeff = F(0,0);
      zeros(0,0) = 1;
    }
  }

  //final = _lambda * (coeff * A + (J - 1) * zeros);
  final = _lambda * (coeff * A + (J - 1 - _ratio) * zeros);
  //final = _lambda * (coeff * A);
  return final;
}
