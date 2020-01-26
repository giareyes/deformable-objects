#ifndef SCALED_NEO_HOOKEAN_H
#define SCALED_NEO_HOOKEAN_H

#include "MATERIAL.h"

class SCALED_NEO_HOOKEAN : public MATERIAL
{
public:
  SCALED_NEO_HOOKEAN(const Real lambda = 4, const Real mu = 1);
  ~SCALED_NEO_HOOKEAN() {};

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

  Real _alpha;
};

#endif
