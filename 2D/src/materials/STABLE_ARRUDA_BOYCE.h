#ifndef STABLE_ARRUDA_BOYCE_H
#define STABLE_ARRUDA_BOYCE_H

#include "MATERIAL.h"

class STABLE_ARRUDA_BOYCE : public MATERIAL
{
public:
  STABLE_ARRUDA_BOYCE(const Real lambda = 1, const Real mu = 1);
  ~STABLE_ARRUDA_BOYCE() {};

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

  Real _beta;
  Real _alpha;
};

#endif
