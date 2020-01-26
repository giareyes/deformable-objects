#ifndef VOLUME_MR_H
#define VOLUME_MR_H

#include "MATERIAL.h"

///////////////////////////////////////////////////////////////////////
// Mooney-Rivlin compression resistence,
// Psi = (lambda / 2) (J - 1)^2
///////////////////////////////////////////////////////////////////////
class VOLUME_MR : public MATERIAL
{
public:
  VOLUME_MR(const Real lambda = 1000, const Real mu = 5000);
  ~VOLUME_MR() {};

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
  MATRIX2 DPvDF(const MATRIX2& F, int i, int j);

  Real _lambda;
  Real _mu;

  Real _ratio;
};

#endif
