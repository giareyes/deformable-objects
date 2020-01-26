#ifndef SYMMETRIC_ARAP_H
#define SYMMETRIC_ARAP_H

#include "MATERIAL.h"

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
class SYMMETRIC_ARAP : public MATERIAL
{
public:
  SYMMETRIC_ARAP(const Real lambda = 1, const Real mu = 1);
  ~SYMMETRIC_ARAP() {};

  // P = first Piola-Kirchoff stress tensor
  // P = F * S
  MATRIX PK1(const MATRIX2& F);

  // S = second Piola-Kirchoff stress tensor
  MATRIX PK2(const MATRIX2& F);
  
  // derivative of PK1 w.r.t. F
  MATRIX DPDF(const MATRIX& F);

  // flatten a matrix into a vector, stacking each of the columns
  // on top of each other
  static VECTOR flatten(const MATRIX& A);

  // get the strain energy
  Real psi(const MATRIX2& F);

  // get the eigenvalues of DPDF, if we have expressions for them
  VEC4 eigenvalues(const MATRIX2& F);

  // get the eigensysem of DPDF, if we have expressions for them
  void eigensystem(const MATRIX2& F, MATRIX4& vectors, VEC4& values);

private:
  // take the svd of F, so we can then form the polar decomposition
  void svd(const MATRIX2& F, MATRIX2& U, VECTOR& S, MATRIX2& V) const;

  // take the polar decomposition, using SVD
  void polarDecomposition(const MATRIX2& F, MATRIX2& R, MATRIX2& S) const;
  
  // take the rotation variant SVD
  void rotationVariantSVD(const MATRIX2& F, MATRIX2& U, MATRIX2& sigma, MATRIX2& V) const;

  Real _lambda;
  Real _mu;
};

#endif
