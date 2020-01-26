#ifndef HYPER_TAN_H
#define HYPER_TAN_H

#include "MATERIAL.h"

///////////////////////////////////////////////////////////////////////
// Hyper-tangent interpolant that ramps from mu0 to mu1
// based on a threshold in the input material Psi
//
// Psi = mu(psi) * psi
//
// mu(psi) = mu0 * (1 - T(psi)) + mu1 * T(psi)
// T = 0.5 + 0.5 * tanh(k * (psi - threshold))
///////////////////////////////////////////////////////////////////////
class HYPER_TAN : public MATERIAL
{
public:
  HYPER_TAN(const Real mu0, const Real mu1, MATERIAL* psi1D);
  ~HYPER_TAN() {};

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

  void setK(const Real k) { _k = k; };
  void setHardeningLength(const Real length) {
    _threshold = pow(length - 1.0, 2.0);
  };

//private:
  MATRIX DmuDF(const MATRIX2& F);
  Real DmuDF(const MATRIX2& F, int index);
  Real mu(const MATRIX2& F) const;

  MATERIAL* _psi1D;
  // sharpness of the tanh transition
  Real _k;

  // transition threshold
  Real _threshold;
  Real _mu0;
  Real _mu1;
};

#endif
