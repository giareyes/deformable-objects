#ifndef STVK_H
#define STVK_H

#include "MATERIAL.h"
#include "EXTRAFUNCTIONS.h"

// finite difference method for DPDF
void HessianDifference();

// finite difference method for PK1
void PK1Difference();

class STVK : public MATERIAL
{
public:
  // settings from the 2008 cubature paper
  STVK(const Real lambda = 1000, const Real mu = 5000);
  ~STVK() {};

  // P = first Piola-Kirchoff stress tensor
  MATRIX PK1(const MATRIX2& F);

  // derivative of PK1 w.r.t. F
  MATRIX DPDF(const MATRIX& F);

  // get the strain energy
  Real psi(const MATRIX2& F);

  Real getLambda();

  Real getMu();

protected:
  Real _lambda;
  Real _mu;
};

#endif
