#include "INVERTIBLE.h"
#include <iostream>

using namespace std;

INVERTIBLE::INVERTIBLE(MATERIAL* material, const Real& eps) :
  _material(material), _eps(eps)
{
  _name = material->name() + string("-INVTBL");
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX INVERTIBLE::PK1(const MATRIX2& F, const VEC2& oldDirection, const bool oldClamped, 
                       VEC2& newDirection, bool& newClamped)
{
  MATRIX2 U, Fhat, FFiltered, V;
  svd(F, oldDirection, oldClamped, U, Fhat, FFiltered, V, newDirection, newClamped);
  MATRIX2 filteredPK1 = _material->PK1(FFiltered);
  return U * filteredPK1 * V.transpose();
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX INVERTIBLE::PK1(const MATRIX2& F)
{
  MATRIX2 U, Fhat, FFiltered, V;
  svd(F, U, Fhat, FFiltered, V);
  MATRIX2 filteredPK1 = _material->PK1(FFiltered);
  return U * filteredPK1 * V.transpose();
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX INVERTIBLE::PK2(const MATRIX2& F, const VEC2& oldDirection, const bool oldClamped,
                       VEC2& newDirection, bool& newClamped)
{
  return F.inverse() * PK1(F, oldDirection, oldClamped, newDirection, newClamped);
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX INVERTIBLE::PK2(const MATRIX2& F)
{
  return F.inverse() * PK1(F);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX INVERTIBLE::DPDF(const MATRIX& F, const VEC2& oldDirection, const bool oldClamped,
                        VEC2& newDirection, bool& newClamped)
{
  MATRIX2 U, Fhat, FFiltered, V;
  svd(F, oldDirection, oldClamped, U, Fhat, FFiltered, V, newDirection, newClamped);

  MATRIX DPFilteredDF = _material->DPDF(FFiltered);

  MATRIX final(4,4);
  int index = 0;
  for (int j = 0; j < 2; j++)
    for (int i = 0; i < 2; i++, index++)
    {
      MATRIX2 DFDF;
      DFDF.setZero();
      DFDF(i,j) = 1;
      MATRIX2 DFhatDF = U.transpose() * DFDF * V;
      VECTOR flatDFhatDF = flatten(DFhatDF);

      VECTOR DPFilteredDFij = DPFilteredDF * flatDFhatDF;

      MATRIX DPDFij = U * unflatten(DPFilteredDFij) * V.transpose();

      final.col(index) = flatten(DPDFij);
    }

  return final;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX INVERTIBLE::DPDF(const MATRIX& F)
{
  MATRIX2 U, Fhat, FFiltered, V;
  svd(F, U, Fhat, FFiltered, V);

  MATRIX DPFilteredDF = _material->DPDF(FFiltered);

  MATRIX final(4,4);
  int index = 0;
  for (int j = 0; j < 2; j++)
    for (int i = 0; i < 2; i++, index++)
    {
      MATRIX2 DFDF;
      DFDF.setZero();
      DFDF(i,j) = 1;
      MATRIX2 DFhatDF = U.transpose() * DFDF * V;
      VECTOR flatDFhatDF = flatten(DFhatDF);

      VECTOR DPFilteredDFij = DPFilteredDF * flatDFhatDF;

      MATRIX DPDFij = U * unflatten(DPFilteredDFij) * V.transpose();

      final.col(index) = flatten(DPDFij);
    }

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
Real INVERTIBLE::psi(const MATRIX2& F)
{
  MATRIX2 U, Fhat, FFiltered, V;
  svd(F, U, Fhat, FFiltered, V);

  return _material->psi(Fhat);
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
Real INVERTIBLE::psi(const MATRIX2& F, const VEC2& oldDirection, const bool oldClamped,
                     VEC2& newDirection, bool& newClamped)
{
  MATRIX2 U, Fhat, FFiltered, V;
  svd(F, oldDirection, oldClamped, U, Fhat, FFiltered, V, newDirection, newClamped);

  return _material->psi(Fhat);
}

///////////////////////////////////////////////////////////////////////
// take the svd of F, so we can form the polar decomposition
///////////////////////////////////////////////////////////////////////
void INVERTIBLE::svd(const MATRIX2& F, MATRIX2& U, MATRIX2& Fhat, MATRIX2& FFiltered, MATRIX2& V)
{
  JacobiSVD<MATRIX> svd(F, ComputeFullU | ComputeFullV);
  U = svd.matrixU();
  V = svd.matrixV();
  VECTOR s = svd.singularValues();
  
  Fhat << s[0], 0,
          0,    s[1];

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
  Fhat = Fhat * L;

  FFiltered = Fhat;
  for (int x = 0; x < 2; x++)
    if (FFiltered(x,x) < _eps)
      FFiltered(x,x) = _eps;
}

///////////////////////////////////////////////////////////////////////
// take the svd of F, so we can form the polar decomposition
///////////////////////////////////////////////////////////////////////
void INVERTIBLE::svd(const MATRIX2& F, const VEC2& oldDirection, const bool oldClamped, MATRIX2& U, MATRIX2& Fhat, MATRIX2& FFiltered, MATRIX2& V, VEC2& newDirection, bool& newClamped)
{
  JacobiSVD<MATRIX> svd(F, ComputeFullU | ComputeFullV);
  U = svd.matrixU();
  V = svd.matrixV();
  VECTOR s = svd.singularValues();
  
  Fhat << s[0], 0,
          0,    s[1];

  // don't use the minimal singular value, use the column of V^T
  // (row of V) that is closest to the previous reflection direction
  VEC2 reflectionDots = V * oldDirection;
  reflectionDots = reflectionDots.array().abs();
  int whichColumn = (reflectionDots[0] > reflectionDots[1]) ? 0 : 1;

  /*
  // if the previous state was not clamped to begin with, ignore
  // the previous reflection direction
  if (oldClamped == false)
    whichColumn = 1;
    */
  whichColumn = 1;

  // load up the reflection
  MATRIX2 L;
  L = MATRIX2::Identity();
  L(whichColumn, whichColumn) = U.determinant() * V.determinant();

  // store the direction we used
  newDirection = V.transpose().col(whichColumn);

  const Real detU = U.determinant();
  const Real detV = V.determinant();

  // which to pull the reflection out of?
  if (detU < 0 && detV > 0)
    U = U * L;
  else if (detU > 0 && detV < 0)
    V = V * L;
  Fhat = Fhat * L;

  FFiltered = Fhat;
  newClamped = false;
  for (int x = 0; x < 2; x++)
    if (FFiltered(x,x) < _eps)
    {
      FFiltered(x,x) = _eps;
      newClamped = true;
    }
}
/*
{
  JacobiSVD<MATRIX> svd(F, ComputeFullU | ComputeFullV);
  U = svd.matrixU();
  V = svd.matrixV();
  VECTOR s = svd.singularValues();
  
  Fhat << s[0], 0,
          0,    s[1];
  
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
  Fhat = Fhat * L;

  FFiltered = Fhat;
  for (int x = 0; x < 2; x++)
    if (FFiltered(x,x) < _eps)
      FFiltered(x,x) = _eps;
}
*/

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR INVERTIBLE::flatten(const MATRIX& A)
{
  VECTOR final(A.rows() * A.cols());

  int i = 0;
  for (int y = 0; y < A.cols(); y++)
    for (int x = 0; x < A.rows(); x++, i++)
      final[i] = A(x,y);

  return final;
}

///////////////////////////////////////////////////////////////////////
// unflatten a vector into a matrix, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
MATRIX2 INVERTIBLE::unflatten(const VECTOR& v)
{
  // make sure the vector size doesn't leave leftovers
  assert(v.size() == 4);
  MATRIX2 final;

  int i = 0;
  for (int y = 0; y < 2; y++)
    for (int x = 0; x < 2; x++, i++)
      final(x,y) = v[i]; 

  return final;
}
