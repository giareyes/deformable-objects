#include "NEO_HOOKEAN.h"
#include <iostream>

using namespace std;

NEO_HOOKEAN::NEO_HOOKEAN(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  _name = string("Neo-Hookean");
  std::cout << " Neo-Hookean settings: " << std::endl;
  std::cout << " \t lambda: " << _lambda << std::endl;
  std::cout << " \t mu:     " << _mu << std::endl;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX NEO_HOOKEAN::PK1(const MATRIX2& F)
{
  //return F * PK2(F);

  const MATRIX2 FTinv = F.transpose().inverse();
  const Real J = F.determinant();

  return _mu * (F - FTinv) + _lambda * log(J) * FTinv;
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX NEO_HOOKEAN::PK2(const MATRIX2& F)
{
  MATRIX2 C = F.transpose() * F;
  MATRIX2 Cinv = C.inverse();
  Real J = F.determinant();

  MATRIX2 S = _mu * (MATRIX2::Identity() - Cinv) + (_lambda * log(J)) * Cinv;
  return S;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX NEO_HOOKEAN::DPDF(const MATRIX& F)
{
  const Real detF = F.determinant();
  const Real invDetF = 1.0 / detF;
  const Real invDetFSq = invDetF * invDetF;

  // inverse transpose of F
  MATRIX2 FinvT;
  FinvT <<  F(1,1), -F(1,0),
           -F(0,1),  F(0,0);
  FinvT *= invDetF;

  // get each partial w.r.t FinvT
  MATRIX2 DFDF[4];
  MATRIX2& DFDF00 = DFDF[0];
  MATRIX2& DFDF10 = DFDF[1];
  MATRIX2& DFDF01 = DFDF[2];
  MATRIX2& DFDF11 = DFDF[3];

  const Real F01F10 = F(0,1) * F(1,0);
  const Real F01F11 = F(0,1) * F(1,1);
  const Real F10F11 = F(1,0) * F(1,1);
  const Real F00F11 = F(0,0) * F(1,1);
  const Real F00F10 = F(0,0) * F(1,0);
  const Real F00F01 = F(0,0) * F(0,1);
  
  DFDF00 << -F(1,1) * F(1,1),  F10F11,
             F01F11,          -F01F10;

  DFDF10 <<  F01F11,          -F00F11,
            -F(0,1) * F(0,1),  F00F01;

  DFDF01 <<  F10F11,          -F(1,0) * F(1,0),
            -F00F11,           F00F10;
  
  DFDF11 << -F01F10,           F00F10,
             F00F01,          -F(0,0) * F(0,0);

  DFDF00 *= invDetFSq;
  DFDF01 *= invDetFSq;
  DFDF10 *= invDetFSq;
  DFDF11 *= invDetFSq;

  const Real lambdaTerm = _lambda * log(detF) - _mu;
  MATRIX final(4,4);

  //const Real logs[] = {log(F(1,1)), log(-F(0,1)), log(-F(1,0)), log(F(0,0))};
  const Real logs[] = {F(1,1) * invDetF, -F(0,1) * invDetF, -F(1,0) * invDetF, F(0,0) * invDetF};

  MATRIX2 zeros;
  zeros.setZero();
  int i = 0;
  for (int y = 0; y < 2; y++)
    for (int x = 0; x < 2; x++, i++)
    {
      zeros(x,y) = _mu;
      MATRIX2 DPDF = zeros + lambdaTerm * DFDF[i] + _lambda * logs[i] * FinvT;
      final.col(i) = flatten(DPDF);
      zeros(x,y) = 0;

      /*
      cout << "=======================================" << endl;
      cout << " xy: " << x << " " << y << endl;
      cout << "=======================================" << endl;
      cout << " DFDF: " << endl;
      cout << DFDF[i] << endl;
      cout << " log term: " << logs[i] << endl;
      */
    }

  return final;
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR NEO_HOOKEAN::flatten(const MATRIX& A)
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
Real NEO_HOOKEAN::psi(const MATRIX2& F)
{
  MATRIX2 C = F.transpose() * F;
  Real Ic = C.trace();
  Real J = F.determinant();
  Real logJ = log(J);

  return _mu * 0.5 * (Ic - 3) - _mu * logJ + _lambda * 0.5 * logJ * logJ;
}
