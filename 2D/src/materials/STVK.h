#ifndef STVK_H
#define STVK_H

#include "MATERIAL.h"

class STVK : public MATERIAL
{
public:
  // settings from the 2008 cubature paper
  STVK(const Real lambda = 1000, const Real mu = 5000);
  ~STVK() {};

  // P = first Piola-Kirchoff stress tensor
  // P = F * S
  MATRIX PK1(const MATRIX2& F);

  // S = second Piola-Kirchoff stress tensor
  MATRIX PK2(const MATRIX2& F);
  
  // derivative of PK2 w.r.t E = 0.5 * (F^T * F - I)
  MATRIX DSDE();

  // derivative of PK1 w.r.t. F
  MATRIX DPDF(const MATRIX& F);

  // get the strain energy
  Real psi(const MATRIX2& F);

protected:
  // take the gradient of FS (i.e. PK1) w.r.t F assuming that S is frozen
  // \frac{\partial F S}{\partial F}
  MATRIX DFSDF(const MATRIX& S);

  // repeat the given matrix "repeats" times along the block diagonal
  MATRIX blockDiag(const MATRIX& A, const int repeats);

  // compute the derivative of E = 1/2 (F^T * F - I)
  // with respect to F
  MATRIX computeDEDF(const MATRIX2& F);

  virtual Real constant() const;
  virtual Real length(const Real sigma) const;
  virtual Real area(const Real s1s2) const;

  Real _lambda;
  Real _mu;
};

#endif
