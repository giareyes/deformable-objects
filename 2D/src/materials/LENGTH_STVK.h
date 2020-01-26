#ifndef LENGTH_STVK_H
#define LENGTH_STVK_H

#include "MATERIAL.h"

///////////////////////////////////////////////////////////////////////
// StVK length resistence,
// Psi = mu E : E
///////////////////////////////////////////////////////////////////////
class LENGTH_STVK : public MATERIAL
{
public:
  LENGTH_STVK(const Real lambda = 1000, const Real mu = 5000);
  ~LENGTH_STVK() {};

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
  MATRIX2 DEDF(const MATRIX2& F, int i, int j);

  Real _lambda;
  Real _mu;
};

#endif
