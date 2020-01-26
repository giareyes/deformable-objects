#include "ARRUDA_BOYCE.h"
#include "VOLUME_NH.h"
#include "VOLUME_MR.h"
#include <iostream>

using namespace std;

ARRUDA_BOYCE::ARRUDA_BOYCE(const Real lambda, const Real mu) :
  _volume(NULL)
{
  _name = string("Arruda-Boyce");

  _volume = new VOLUME_NH(lambda, mu);
  //_volume = new VOLUME_MR(lambda, mu);

  //_N = 5;
  //_N = 2;
  //_N = 4;
  _N = 1;
  //_N = 0.1;
  //_N = 100;
  //_C1 = 5000;
  _C1 = mu;
  //_K = 100000;
  _K = lambda;

  cout << " Arruda-Boyce settings, N = "  << _N << " C1 = " << _C1 << " K = " << _K << endl;
  //exit(0);
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX ARRUDA_BOYCE::PK1(const MATRIX2& F)
{
  const MATRIX2 C = F.transpose() * F;
  const Real Ic = C.trace();
  const Real Ic2 = Ic * Ic;
  const Real Ic3 = Ic2 * Ic;
  const Real Ic4 = Ic3 * Ic;

  const Real N2 = _N * _N;
  const Real N3 = N2 * _N;
  const Real N4 = N3 * _N;

  const Real coeff = _C1 * (1.0 + 
                            1.0 / (5.0 * _N) * Ic + 
                            11.0 / (175.0 * N2) * Ic2 +
                            19.0 / (875.0 * N3) * Ic3 +
                            519.0 / (67375 * N4) * Ic4);

  if (_volume)
    return coeff * F + _volume->PK1(F);
  else
    return coeff * F;
  /*
  //return coeff * F;

  // add a volume penalty term
  Real J = F.determinant();
  MATRIX2 DJDF;
  DJDF <<  F(1,1), -F(1,0),
          -F(0,1), F(0,0);

  return coeff * F + _K * (J - 1) * DJDF;
  */
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX ARRUDA_BOYCE::PK2(const MATRIX2& F)
{
  return F.inverse() * PK1(F);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX ARRUDA_BOYCE::DPDF(const MATRIX& F)
{
  const MATRIX2 C = F.transpose() * F;
  const Real Ic = C.trace();
  const Real Ic2 = Ic * Ic;
  const Real Ic3 = Ic2 * Ic;
  const Real Ic4 = Ic3 * Ic;

  const Real N2 = _N * _N;
  const Real N3 = N2 * _N;
  const Real N4 = N3 * _N;

  const Real coeff = _C1 * (1.0 + 
                            1.0 / (5.0 * _N) * Ic + 
                            11.0 / (175.0 * N2) * Ic2 +
                            19.0 / (875.0 * N3) * Ic3 +
                            519.0 / (67375 * N4) * Ic4);
  const Real dCoeff = _C1 * (2.0 / (5.0 * _N) + 
                             44.0 / (175.0 * N2) * Ic +
                             114.0 / (875.0 * N3) * Ic2 +
                             4152.0 / (67375 * N4) * Ic3);

  MATRIX final(4,4);
  int index = 0; 
  for (int j = 0; j < 2; j++)
    for (int i = 0; i < 2; i++, index++)
    {
      const MATRIX2 DFDFij = DFDF(F, i, j);
      MATRIX2 column = coeff * DFDFij + dCoeff * F(i,j) * F;
      //column += DPvDF(F, i, j);

      final.col(index) = flatten(column);
    }

  if (_volume)
    final += _volume->DPDF(F);

  return final;
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR ARRUDA_BOYCE::flatten(const MATRIX& A)
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
Real ARRUDA_BOYCE::psi(const MATRIX2& F)
{
  /*
  MATRIX2 C = F.transpose() * F;
  Real Ic = C.trace();
  Real IIc = (C * C).trace();

  Real J = F.determinant();

  Real volume = _lambda * 0.5 * (J - 1) * (J - 1);

  //cout << " Ic: " << Ic << endl;

  //return _mu10 * (Ic - 3) + _mu01 * 0.5 * (Ic * Ic - IIc - 3);
  //return _mu10 * (Ic - 2) + _mu01 * 0.5 * (Ic * Ic - IIc - 2);
  return _mu10 * (Ic - 2) + _mu01 * 0.5 * (Ic * Ic - IIc - 2) + volume;
  */
  MATRIX2 C = F.transpose() * F;
  Real Ic = C.trace();
  
  //Real J = F.determinant();
  //Real volume = _K * 0.5 * (J - 1) * (J - 1);
  Real volume = 0;

  if (_volume)
    volume += _volume->psi(F);

  return _C1 * (0.5 * (Ic - 3) +
                1.0 / (20.0 * _N) *              (Ic * Ic - 9) +
                11.0 / (1050.0 * _N * _N) *      (Ic * Ic * Ic - 27) +
                19.0 / (7000.0 * _N * _N * _N) * (pow(Ic, 4.0) - 81) +
                519.0 / (673750 * pow(_N, 4)) *  (pow(Ic, 5.0) - 243)) + volume;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX2 ARRUDA_BOYCE::DFDF(const MATRIX2& F, int i, int j)
{
  MATRIX2 zero;
  zero.setZero();
  zero(i,j) = 1;
  return zero;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX2 ARRUDA_BOYCE::DPvDF(const MATRIX2& F, int i, int j)
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

  final = _K * (coeff * A + (J - 1) * zeros);
  return final;
}
