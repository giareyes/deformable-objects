#ifndef COROTATIONAL_H
#define COROTATIONAL_H

#include "MATERIAL.h"

class COROTATIONAL : public MATERIAL
{
public:
  COROTATIONAL(const Real lambda = 1000, const Real mu = 5000);
  ~COROTATIONAL() {};

  // P = first Piola-Kirchoff stress tensor
  // P = F * S
  MATRIX PK1(const MATRIX2& F);

  // S = second Piola-Kirchoff stress tensor
  MATRIX PK2(const MATRIX2& F);
  
  // derivative of PK1 w.r.t. F
  MATRIX DPDF(const MATRIX& F);

  // derivative of rigid rotation w.r.t. F
  MATRIX DRDF(const MATRIX& F);
  
  // derivative of rigid scaling w.r.t. F
  MATRIX DSDF(const MATRIX& F);

  // nearest rotation
  MATRIX2 R(const MATRIX& F);
  
  // stretching-only term
  MATRIX2 S(const MATRIX& F);

  // take the svd of F, so we can then form the polar decomposition
  void svd(const MATRIX2& F, MATRIX2& U, VECTOR& S, MATRIX2& V) const;
  void svd(const MATRIX2& F, MATRIX2& U, MATRIX2& S, MATRIX2& V) const;

  // get the gradient of each component of the SVD
  void svdGradient(const MATRIX2& F, const int i, const int j,
                   MATRIX2& DU, MATRIX2& DS, MATRIX2& DV);

  // run a unit test on the polar gradients
  void testPolarGrad(const MATRIX2& F);

  // get the strain energy
  Real psi(const MATRIX2& F);

private:
  // take the polar decomposition, using SVD
  void polarDecomposition(const MATRIX2& F, MATRIX2& R, MATRIX2& S) const;

  // get the gradient of the scaling
  MATRIX DSDF(const MATRIX2& F, int index);
  MATRIX DSDF(const MATRIX2& F, const MATRIX2& R, const MATRIX2& S, const MATRIX2& DR, int index);
  
  // get the gradient of the rotation
  MATRIX DRDF(const MATRIX2& F, int index);
  MATRIX DRDF(const MATRIX2& F, const MATRIX2& R, const MATRIX2& S, int index);
  MATRIX DRTDF(const MATRIX2& F, int index);

  // get the gradient of F w.r.t PK2, got a specific entry of F
  MATRIX2 DPDF(const MATRIX2& F, int index);
  MATRIX2 DPDFFast(const MATRIX2& F, int index);

  // promote a 2x2 matrix to 3x3, with zero off-diagonals and a
  // and along the diagonal
  MATRIX3 promote(const MATRIX2& A);
  
  // demote a 3x3 to a 2x2, take the upper left of the matrix
  MATRIX2 demote(const MATRIX3& A);

  // get the skew-symmetric portion of this matrix
  VECTOR skew(const MATRIX3& A);

  // flatten a matrix into a vector, stacking each of the columns
  // on top of each other
  static VECTOR flatten(const MATRIX& A);

  Real _lambda;
  Real _mu;
};

#endif
