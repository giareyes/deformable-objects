#ifndef MOONEY_RIVLIN_H
#define MOONEY_RIVLIN_H

#include "MATERIAL.h"

class MOONEY_RIVLIN : public MATERIAL
{
public:
  MOONEY_RIVLIN(const Real lambda = 1000, const Real mu = 5000);
  ~MOONEY_RIVLIN() {};

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
  MATRIX2 DFDF(const MATRIX2& F, int i, int j);
  MATRIX2 DCDF(const MATRIX2& F, int i, int j);
  MATRIX2 DPvDF(const MATRIX2& F, int i, int j);

  // the volume preservation model
  MATERIAL* _volume;

  // the passed in parameters
  Real _lambda;
  Real _mu;

  // the actual parameters used, based on what was passed in
  Real _mu01;
  Real _mu10;
};

#endif
