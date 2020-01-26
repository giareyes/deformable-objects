#ifndef INVERTIBLE_H
#define INVERTIBLE_H

#include "MATERIAL.h"

///////////////////////////////////////////////////////////////////////
// invertible material, filters the incoming F to preclude
// badly-behaved regions
///////////////////////////////////////////////////////////////////////
class INVERTIBLE : public MATERIAL
{
public:
  INVERTIBLE(MATERIAL* material, const Real& eps);
  ~INVERTIBLE() {};

  // P = first Piola-Kirchoff stress tensor
  // P = F * S
  MATRIX PK1(const MATRIX2& F);
  MATRIX PK1(const MATRIX2& F, const VEC2& oldDirection, const bool oldClamped, VEC2& newDirection, bool& newClamped);

  // S = second Piola-Kirchoff stress tensor
  MATRIX PK2(const MATRIX2& F);
  MATRIX PK2(const MATRIX2& F, const VEC2& oldDirection, const bool oldClamped, VEC2& newDirection, bool& newClamped);
  
  // derivative of PK1 w.r.t. F
  MATRIX DPDF(const MATRIX& F);
  MATRIX DPDF(const MATRIX& F, const VEC2& oldDirection, const bool oldClamped, VEC2& newDirection, bool& newClamped);

  // get the strain energy
  Real psi(const MATRIX2& F);
  Real psi(const MATRIX2& F, const VEC2& oldDirection, const bool oldClamped, VEC2& newDirection, bool& newClamped);

protected:
  void svd(const MATRIX2& F, MATRIX2& U, MATRIX2& Fhat, MATRIX2& FFiltered, MATRIX2& V);
  void svd(const MATRIX2& F, const VEC2& oldDirection, const bool oldClamped, MATRIX2& U, MATRIX2& Fhat, MATRIX2& FFiltered, MATRIX2& V, VEC2& newDirection, bool& newClamped);
  
  // flatten a matrix into a vector, stacking each of the columns
  // on top of each other
  static VECTOR flatten(const MATRIX& A);

  // unflatten a vector into a matrix, stacking each of the columns
  // on top of each other
  static MATRIX2 unflatten(const VECTOR& v);

  MATERIAL* _material;
  Real _eps;
};

#endif
