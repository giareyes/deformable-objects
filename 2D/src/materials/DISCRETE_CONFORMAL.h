#ifndef DISCRETE_CONFORMAL_H
#define DISCRETE_CONFORMAL_H

#include "MATERIAL.h"

class DISCRETE_CONFORMAL : public MATERIAL
{
public:
  DISCRETE_CONFORMAL(const Real lambda = 20, const Real mu = 1);
  ~DISCRETE_CONFORMAL() {};

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

private:
  // take the polar decomposition, using SVD
  void polarDecomposition(const MATRIX2& F, MATRIX2& R, MATRIX2& S) const;
  
  // take the polar decomposition, using SVD
  void svd(const MATRIX2& F, MATRIX2& U, MATRIX2& sigma, MATRIX2& V) const;

  // clamp eigenvalues to zero
  MATRIX clampEigenvalues(const MATRIX& A);

  Real _lambda;
  Real _mu;

  Real _alpha;
  Real _beta;
};

#endif
