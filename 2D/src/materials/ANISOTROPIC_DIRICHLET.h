#ifndef ANISOTROPIC_DIRICHLET_H
#define ANISOTROPIC_DIRICHLET_H

#include "MATERIAL.h"

///////////////////////////////////////////////////////////////////////
// anisotropic shrink-to-a-point energy
///////////////////////////////////////////////////////////////////////
class ANISOTROPIC_DIRICHLET : public MATERIAL
{
public:
  ANISOTROPIC_DIRICHLET(const Real lambda = 1000, const Real mu = 5000);
  ~ANISOTROPIC_DIRICHLET() {};

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

  static bool testSVD();

  void setAnisotropyDirection(const VEC2& a);

  const VEC2 anisotropyDirection() const { return _a; };

private:
  // promote a 2x2 matrix to 3x3, with zero off-diagonals and a
  // and along the diagonal
  MATRIX3 promote(const MATRIX2& A);

  // demote a 3x3 to a 2x2, take the upper left of the matrix
  MATRIX2 demote(const MATRIX3& A);

  // get the skew-symmetric portion of this matrix
  VECTOR skew(const MATRIX3& A);

  // take the svd of F, so we can then form the polar decomposition
  void mySVD(const MATRIX2& F, MATRIX2& U, VECTOR& S, MATRIX2& V) const;
  void svd(const MATRIX2& F, MATRIX2& U, VECTOR& S, MATRIX2& V) const;

  // take the polar decomposition, using SVD
  void polarDecomposition(const MATRIX2& F, MATRIX2& R, MATRIX2& S) const;

  // get the gradient of the scaling
  MATRIX DSDF(const MATRIX2& F, const MATRIX2& R, const MATRIX2& S, const MATRIX2& DR, int index);

  // get the gradient of the rotation
  MATRIX DRDF(const MATRIX2& F, const MATRIX2& R, const MATRIX2& S, int index);

  // generate coefficients for Jacobi (Givens) rotation
  bool symSchur2(const MATRIX2& A, Real& c, Real& s) const;

  Real _lambda;
  Real _mu;

  // anisotropy direction
  VEC2 _a;
  MATRIX2 _A;
};

#endif
