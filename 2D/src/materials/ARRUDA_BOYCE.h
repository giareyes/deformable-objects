#ifndef ARRUDA_BOYCE_H
#define ARRUDA_BOYCE_H

#include "MATERIAL.h"

class ARRUDA_BOYCE : public MATERIAL
{
public:
  ARRUDA_BOYCE(const Real lambda = 1000, const Real mu = 5000);
  ~ARRUDA_BOYCE() { if (_volume) delete _volume; };

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
  MATRIX2 DPvDF(const MATRIX2& F, int i, int j);

  // material parameters
  Real _C1;
  Real _N;
  Real _K;
  
  // volume preservation term
  MATERIAL* _volume;
};

#endif
