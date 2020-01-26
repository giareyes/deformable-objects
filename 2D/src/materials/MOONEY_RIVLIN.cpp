#include "MOONEY_RIVLIN.h"
#include "VOLUME_MR.h"
#include "VOLUME_NH.h"
#include <iostream>

using namespace std;

MOONEY_RIVLIN::MOONEY_RIVLIN(const Real lambda, const Real mu) :
  _volume(NULL), _lambda(lambda), _mu(mu)
{
  _name = string("Mooney-Rivlin");

  _volume = new VOLUME_MR(_lambda, _mu);
  //_volume = new VOLUME_NH(_lambda, _mu);

  // convert back to E and nu
  const Real E = _mu * (3.0 * _lambda + 2.0 * _mu) / (_lambda + _mu);
  const Real nu = _lambda / (2.0 * (_lambda + _mu));

  const Real G = E / (2.0 * (1 + nu));
  //_mu10 = G / (2.0 * (1.0 + mu));
  //_mu01 = mu * _mu10;
  
  _mu10 = _mu * 0.5;
  //_mu01 = mu * 2;
  _mu01 = _mu10 * 0.5;
  
  //_mu10 = mu;
  //_mu01 = 0;
 
  // some measured values found online 
  //_mu10 = 0.035e6;
  //_mu01 = 0.013e6;
  
  cout << " Mooney-Rivlin settings mu10: "  << _mu10 << " mu01: " << _mu01 << " lambda: " << _lambda << endl;
  //exit(0);
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX MOONEY_RIVLIN::PK1(const MATRIX2& F)
{
  MATRIX2 C = F.transpose() * F;
  Real Ic = C.trace();
#if 0
  Real J = F.determinant();

  MATRIX2 DJDF;
  DJDF <<  F(1,1), -F(1,0),
          -F(0,1), F(0,0);

  return 2.0 * F * ((_mu10 + _mu01 * Ic) * MATRIX2::Identity() - _mu01 * C) + _lambda * (J - 1) * DJDF; 
#else
  //return 2.0 * F * ((_mu10 + _mu01 * Ic) * MATRIX2::Identity() - _mu01 * C) + _volume->PK1(F);
  return 2.0 * (_mu10 + _mu01 * Ic) * F - 2.0 * _mu01 * F * C + _volume->PK1(F);
#endif
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX MOONEY_RIVLIN::PK2(const MATRIX2& F)
{
  /*
  MATRIX2 C = F.transpose() * F;
  Real Ic = C.trace();

  MATRIX2 S = 2.0 * ((_mu10 + _mu01 * Ic) * MATRIX2::Identity() - _mu01 * C); 
  return S;
  */
  return F.inverse() * PK1(F);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX MOONEY_RIVLIN::DPDF(const MATRIX& F)
{
  MATRIX final(4,4);
  const MATRIX2 C = F.transpose() * F;
  const MATRIX2 I = MATRIX2::Identity();
  const MATRIX2 IcI = C.trace() * I;

#if 0
  int index = 0; 
  for (int j = 0; j < 2; j++)
    for (int i = 0; i < 2; i++, index++)
    {
      const MATRIX2 DFDFij = DFDF(F, i, j);
      MATRIX2 column = 2.0 * (_mu10 * DFDFij +
                              _mu01 * DFDFij * (IcI  - C) +
                              _mu01 * F * (2 * F(i,j) * I - DCDF(F,i,j)));

      column += DPvDF(F, i, j);
      final.col(index) = flatten(column);
    }

  //return final + _volume->DPDF(F);
  return final;
#else
  int index = 0; 
  for (int j = 0; j < 2; j++)
    for (int i = 0; i < 2; i++, index++)
    {
      const MATRIX2 DFDFij = DFDF(F, i, j);
      MATRIX2 column = 2.0 * (_mu10 * DFDFij +
                              _mu01 * DFDFij * (IcI  - C) +
                              _mu01 * F * (2 * F(i,j) * I - DCDF(F,i,j)));

      final.col(index) = flatten(column);
    }
  return final + _volume->DPDF(F);
#endif
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR MOONEY_RIVLIN::flatten(const MATRIX& A)
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
Real MOONEY_RIVLIN::psi(const MATRIX2& F)
{
  MATRIX2 C = F.transpose() * F;
  Real Ic = C.trace();
  Real IIc = (C * C).trace();

#if 0
  Real J = F.determinant();
  Real volume = _lambda * 0.5 * (J - 1) * (J - 1);

  //cout << " Ic: " << Ic << endl;

  //return _mu10 * (Ic - 3) + _mu01 * 0.5 * (Ic * Ic - IIc - 3);
  //return _mu10 * (Ic - 2) + _mu01 * 0.5 * (Ic * Ic - IIc - 2);
  return _mu10 * (Ic - 2) + _mu01 * 0.5 * (Ic * Ic - IIc - 2) + volume;
#else
  //return _mu10 * (Ic - 2) + _mu01 * 0.5 * (Ic * Ic - IIc - 2) + _volume->psi(F);
  return _mu10 * (Ic - 3) + _mu01 * 0.5 * (Ic * Ic - IIc - 3) + _volume->psi(F);
#endif
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX2 MOONEY_RIVLIN::DFDF(const MATRIX2& F, int i, int j)
{
  MATRIX2 zero;
  zero.setZero();
  zero(i,j) = 1;
  return zero;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX2 MOONEY_RIVLIN::DPvDF(const MATRIX2& F, int i, int j)
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

  final = _lambda * (coeff * A + (J - 1) * zeros);
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX2 MOONEY_RIVLIN::DCDF(const MATRIX2& F, int i, int j)
{
  MATRIX2 zero;

  if (i == 0)
  {
    // (0,0)
    if (j == 0)
      zero << 2.0 * F(0,0), F(0,1),
                    F(0,1), 0;
    // (0,1)
    else
      zero << 0, F(0,0),
              F(0,0), 2.0 * F(0,1);
  }
  else
  {
    // (1,0)
    if (j == 0)
      zero << 2.0 * F(1,0), F(1,1),
                    F(1,1), 0;
    // (1,1)
    else
      zero << 0, F(1,0),
              F(1,0), 2.0 * F(1,1);
  }
  return zero;
}
