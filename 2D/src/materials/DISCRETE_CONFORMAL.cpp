#include "DISCRETE_CONFORMAL.h"
#include <iostream>

using namespace std;

#define USING_ORIGINAL 1
#define USING_2J 0
#define USING_FUNG 0

DISCRETE_CONFORMAL::DISCRETE_CONFORMAL(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  //_lambda = 1.0;
  //_mu = 1.0;
  //_mu = 0.1;

  _name = string("Discrete_Conformal");
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX DISCRETE_CONFORMAL::PK1(const MATRIX2& F)
{
  MATRIX2 DJ;
  DJ <<  F(1,1), -F(1,0),
        -F(0,1),  F(0,0);
 
#if USING_ORIGINAL 
  MATRIX2 finalReturn = 2 * F - DJ;
#endif
#if USING_2J
  MATRIX2 finalReturn = 2 * F - 2 * DJ;
#endif
#if USING_FUNG
  MATRIX2 R,S;
  polarDecomposition(F,R,S);
  MATRIX2 finalReturn = 2 * F - exp(S.trace() - 2) * R;

  finalReturn *= _mu;
#endif
  return finalReturn;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX DISCRETE_CONFORMAL::PK2(const MATRIX2& F)
{
  return F.inverse() * PK1(F);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX DISCRETE_CONFORMAL::DPDF(const MATRIX& F)
{
  MATRIX finalMatrix(4,4);
  finalMatrix.setIdentity();
  finalMatrix *= 2.0;

  MATRIX anti(4,4);
  anti.setZero();
  anti(0,3) = 1;
  anti(1,2) = -1;
  anti(2,1) = -1;
  anti(3,0) = 1;

#if USING_ORIGINAL 
  finalMatrix -= anti;
#endif
#if USING_2J
  finalMatrix -= 2.0 * anti;
#endif
#if USING_FUNG
  MATRIX2 U,sigma,V;
  svd(F,U,sigma,V);
  MATRIX2 twist;
  twist << 0, -1, 1 ,0;
  Real trS = sigma.trace();

  MATRIX2 Q = (1.0 / sqrt(2.0)) * U * twist * V.transpose();
  VECTOR q = flatten(Q);
  Real eigenvalue = 2.0 / (trS);
  finalMatrix -= exp(trS - 2.0) * eigenvalue * (q * q.transpose());

  MATRIX2 R = U * V.transpose();
  VECTOR r = flatten(R);
  finalMatrix -= exp(trS - 2.0) * (r * r.transpose());
  
  finalMatrix *= _mu;
  //return clampEigenvalues(finalMatrix);
#endif

  return finalMatrix;
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR DISCRETE_CONFORMAL::flatten(const MATRIX& A)
{
  VECTOR final(A.rows() * A.cols());

  int i = 0;
  for (int y = 0; y < A.cols(); y++)
    for (int x = 0; x < A.rows(); x++, i++)
      final[i] = A(x,y);

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the polar decomposition, using SVD
///////////////////////////////////////////////////////////////////////
void DISCRETE_CONFORMAL::svd(const MATRIX2& F, MATRIX2& U, MATRIX2& sigma, MATRIX2& V) const
{
  JacobiSVD<MATRIX> svd(F, ComputeFullU | ComputeFullV);
  U = svd.matrixU();
  V = svd.matrixV();
  VECTOR s = svd.singularValues();

  // pipe the singular values to a matrix
  sigma << s[0], 0,
           0,    s[1];

  MATRIX2 Uoriginal = U;
  MATRIX2 sigmaOriginal = sigma;

  // reflection matrix
  MATRIX2 L;
  L.setZero();
  L(0,0) = 1;
  L(1,1) = (U * V.transpose()).determinant();

  const Real detU = U.determinant();
  const Real detV = V.determinant();

  // which to pull the reflection out of?
  if (detU < 0 && detV > 0)
    U = U * L;
  else if (detU > 0 && detV < 0)
    V = V * L;
  sigma = sigma * L;
}

///////////////////////////////////////////////////////////////////////
// take the polar decomposition, using SVD
///////////////////////////////////////////////////////////////////////
void DISCRETE_CONFORMAL::polarDecomposition(const MATRIX2& F, MATRIX2& R, MATRIX2& S) const
{
  JacobiSVD<MATRIX> svd(F, ComputeFullU | ComputeFullV);
  MATRIX2 U = svd.matrixU();
  MATRIX2 V = svd.matrixV();
  VECTOR s = svd.singularValues();

  // pipe the singular values to a matrix
  MATRIX2 sigma;
  sigma << s[0], 0,
           0,    s[1];

  MATRIX2 Uoriginal = U;
  MATRIX2 sigmaOriginal = sigma;

  // reflection matrix
  MATRIX2 L;
  L.setZero();
  L(0,0) = 1;
  L(1,1) = (U * V.transpose()).determinant();

  const Real detU = U.determinant();
  const Real detV = V.determinant();

  // which to pull the reflection out of?
  if (detU < 0 && detV > 0)
    U = U * L;
  else if (detU > 0 && detV < 0)
    V = V * L;
  sigma = sigma * L;

  R = U * V.transpose();
  S = V * sigma * V.transpose();
}

///////////////////////////////////////////////////////////////////////
// clamp eigenvalues to zero
///////////////////////////////////////////////////////////////////////
MATRIX DISCRETE_CONFORMAL::clampEigenvalues(const MATRIX& A)
{
  EigenSolver<MATRIX> eigensolver(A);
  VECTOR realEigenvalues = eigensolver.eigenvalues().real();
  MATRIX realEigenvectors = eigensolver.eigenvectors().real();

  MATRIX Delta(4,4);
  Delta.setZero();
  Real eps = 1e-4;
  for (int x = 0; x < 4; x++)
    Delta(x,x) = (realEigenvalues[x] > eps) ? realEigenvalues[x] : eps;

  MATRIX clamped = realEigenvectors * Delta * realEigenvectors.inverse();
  return clamped;
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
Real DISCRETE_CONFORMAL::psi(const MATRIX2& F)
{
  Real IC = (F.transpose() * F).trace();
  Real J = F.determinant();

#if USING_ORIGINAL 
  return IC - J;
#endif
#if USING_2J
  return IC - 2 * J;
#endif
#if USING_FUNG
  MATRIX2 R,S;
  polarDecomposition(F,R,S);
  return _mu * (IC - exp(S.trace() - 2.0));
#endif
}
