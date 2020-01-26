#ifndef STABLE_NEO_HOOKEAN_H
#define STABLE_NEO_HOOKEAN_H

#include "MATERIAL.h"

class STABLE_NEO_HOOKEAN : public MATERIAL
{
public:
  STABLE_NEO_HOOKEAN(const Real lambda = 1, const Real mu = 1);
  ~STABLE_NEO_HOOKEAN() {};

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
  Real _lambda;
  Real _mu;
};

#endif
