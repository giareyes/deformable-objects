#include "COROTATIONAL.h"
#include <iostream>
#include "TIMER.h"

using namespace std;

COROTATIONAL::COROTATIONAL(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  _name = string("Co-rotational");
  /*
  std::cout << " Corotational settings: " << std::endl;
  std::cout << " \t lambda: " << _lambda << std::endl;
  std::cout << " \t mu:     " << _mu << std::endl;
  */
}

#define FREEZE_R 0

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX COROTATIONAL::PK1(const MATRIX2& F)
{
  MATRIX2 R,S;
  polarDecomposition(F,R,S);

  MATRIX final;
  const MATRIX2 I = MATRIX2::Identity();
  final = R * (2.0 * _mu * (S - I) + _lambda * (S - I).trace() * I);
  return final;
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX COROTATIONAL::PK2(const MATRIX2& F)
{
  return F.inverse() * PK1(F);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX COROTATIONAL::DPDF(const MATRIX& F)
{
  TIMER functionTimer(__FUNCTION__);
#if FREEZE_R
  MATRIX final(4,4);
  final.setZero();
  for (int i = 0; i < 4; i++)
    final(i,i) = 2.0 * _mu;
#else
  MATRIX final(4,4);
  for (int i = 0; i < 4; i++)
  {
    //MATRIX2 DPDF_i = DPDF(F, i);
    MATRIX2 DPDF_i = DPDFFast(F, i);
    final.col(i) = flatten(DPDF_i);
  }
#endif

  SelfAdjointEigenSolver<MATRIX> eigensolver(final);
  VECTOR numerical = eigensolver.eigenvalues();
  MATRIX eigenvectors = eigensolver.eigenvectors();
  for (int x = 0; x < numerical.size(); x++)
    numerical[x] = (numerical[x] < 0.0) ? 0.0 : numerical[x];
  final = eigenvectors * numerical.asDiagonal() * eigenvectors.transpose();

  /*
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " DRDF: " << endl << DRDF(F) << endl;
  cout << " DSDF: " << endl << DSDF(F) << endl;
  cout << " DPDF: " << endl << final << endl;

  EigenSolver<MATRIX> eigensolver(final);
  //VECTOR values = eigensolver.eigenvalues();
  //MATRIX vectors = eigensolver.eigenvectors();
  cout << " eigenvalues: " << endl << eigensolver.eigenvalues().transpose() << endl;
  //cout << " eigenvectors: " << endl << vectors << endl;
  */

  /*
  MATRIX lambda(4,4);
  lambda.setZero();
  for (int i = 0; i < 4; i++)
  {
    Real value = eigensolver.eigenvalues()[i];
    //lambda(i,i) = (value < 0) ? 0 : value;
    lambda(i,i) = (value > 0) ? -1e-3 : value;
    if (value > 0)
      cout << " CLAMPED " << endl;
  }

  final = vectors * lambda * vectors.inverse();
  */

  return final;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX2 COROTATIONAL::DPDFFast(const MATRIX2& F, int index)
{
  // take the polar decomposition, using SVD
  MATRIX2 R, S;
  polarDecomposition(F,R,S);

  /*
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " F: " << endl << F << endl;
  cout << " index: " << index << endl;
  cout << " R: " << endl << R << endl;
  cout << " S: " << endl << S << endl;
  */

  // lots of optimizations possible here, the polar decomposition
  // and all the gradients have redundant computation
  //MATRIX2 DRDF_i = DRDF(F, index);
  MATRIX2 DRDF_i = DRDF(F, R, S, index);
  MATRIX2 DSDF_i = DSDF(F, R, S, DRDF_i, index);

  /*
  cout << " DRDF: " << endl << DRDF_i << endl;
  cout << " DSDF: " << endl << DSDF_i << endl;
  EigenSolver<MATRIX2> eigensolver(DRDF_i);
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  */

  //DRDF_i << 0,0,0,0;
  //DSDF_i << 0,0,0,1;

  // the necessity for this is a little mysterious. The finite difference
  // verification for DRDF does not bear out the negative, but it is needed
  // for the DPDF and K verifications.
  //
  // SKEW HAS A NEGATIVE BUG IN IT
  DRDF_i *= -1;

  const MATRIX2 I = MATRIX2::Identity();
  MATRIX2 final;
  final = DRDF_i * (2.0 * _mu * (S - I) + _lambda * (S - I).trace() * I);
  final += R * (2.0 * _mu * DSDF_i + _lambda * DSDF_i.trace() * I);
 
  //final.setZero(); 
  //final += DRDF_i * 2.0 * _mu * (S - I);
  //final += R * 2.0 * _mu * DSDF_i;

  //cout << " DSDF: " << endl;
  //cout << DSDF_i << endl;

  return final;
}

#if 0
///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX COROTATIONAL::DPDF(const MATRIX& F)
{
  TIMER functionTimer(__FUNCTION__);
  MATRIX final(4,4);
  for (int i = 0; i < 4; i++)
  {
    //MATRIX2 DPDF_i = DPDF(F, i);
    MATRIX2 DPDF_i = DPDFFast(F, i);
    final.col(i) = flatten(DPDF_i);
  }

  /*
  MATRIX lambda(4,4);
  lambda.setZero();
  for (int i = 0; i < 4; i++)
  {
    Real value = eigensolver.eigenvalues()[i];
    //lambda(i,i) = (value < 0) ? 0 : value;
    lambda(i,i) = (value > 0) ? -1e-3 : value;
    if (value > 0)
      cout << " CLAMPED " << endl;
  }

  final = vectors * lambda * vectors.inverse();
  */

  return final;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX2 COROTATIONAL::DPDFFast(const MATRIX2& F, int index)
{
  // take the polar decomposition, using SVD
  MATRIX2 R, S;
  polarDecomposition(F,R,S);

  // lots of optimizations possible here, the polar decomposition
  // and all the gradients have redundant computation
  MATRIX2 DRDF_i = DRDF(F, R, S, index);
  MATRIX2 DSDF_i = DSDF(F, R, S, DRDF_i, index);

  // the necessity for this is a little mysterious. The finite difference
  // verification for DRDF does not bear out the negative, but it is needed
  // for the DPDF and K verifications.
  DRDF_i *= -1;

  const MATRIX2 I = MATRIX2::Identity();
  MATRIX2 final;
  final = DRDF_i * (2.0 * _mu * (S - I) + _lambda * (S - I).trace() * I);
  final += R * (2.0 * _mu * DSDF_i + _lambda * DSDF_i.trace() * I);

  return final;
}
#endif

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX2 COROTATIONAL::DPDF(const MATRIX2& F, int index)
{
  // lots of optimizations possible here, the polar decomposition
  // and all the gradients have redundant computation
  MATRIX2 DRDF_i = DRDF(F, index);
  MATRIX2 DSDF_i = DSDF(F, index);

  // the necessity for this is a little mysterious. The finite difference
  // verification for DRDF does not bear out the negative, but it is needed
  // for the DPDF and K verifications.
  DRDF_i *= -1;

  // take the polar decomposition, using SVD
  MATRIX2 R, S;
  polarDecomposition(F,R,S);

  const MATRIX2 I = MATRIX2::Identity();
  MATRIX2 final;
  final = DRDF_i * (2.0 * _mu * (S - I) + _lambda * (S - I).trace() * I);
  final += R * (2.0 * _mu * DSDF_i + _lambda * DSDF_i.trace() * I);

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the svd of F, so we can form the polar decomposition
///////////////////////////////////////////////////////////////////////
void COROTATIONAL::svd(const MATRIX2& F, MATRIX2& U, MATRIX2& S, MATRIX2& V) const
{
  JacobiSVD<MATRIX> svd(F, ComputeFullU | ComputeFullV);
  U = svd.matrixU();
  V = svd.matrixV();
  VECTOR Svec = svd.singularValues();

  S.setZero();
  S(0,0) = Svec[0];
  S(1,1) = Svec[1];
}

///////////////////////////////////////////////////////////////////////
// take the svd of F, so we can form the polar decomposition
///////////////////////////////////////////////////////////////////////
void COROTATIONAL::svd(const MATRIX2& F, MATRIX2& U, VECTOR& S, MATRIX2& V) const
{
  JacobiSVD<MATRIX> svd(F, ComputeFullU | ComputeFullV);
  U = svd.matrixU();
  V = svd.matrixV();
  S = svd.singularValues();
}

///////////////////////////////////////////////////////////////////////
// take the polar decomposition, using SVD
///////////////////////////////////////////////////////////////////////
void COROTATIONAL::polarDecomposition(const MATRIX2& F, MATRIX2& R, MATRIX2& S) const
{
  MATRIX2 U;
  VECTOR s;
  MATRIX2 V;
  svd(F, U,s,V);

  // pipe the singular values to a matrix
  MATRIX2 sigma;
  sigma << s[0], 0,
           0,    s[1];

  // reflection matrix
  MATRIX2 L;
  L.setZero();
  L(0,0) = 1;
  L(1,1) = (U * V.transpose()).determinant();
  U = U * L;
  sigma = sigma * L;

  R = U * V.transpose();
  S = V * sigma * V.transpose();
}

///////////////////////////////////////////////////////////////////////
// get the gradient of the rotation
///////////////////////////////////////////////////////////////////////
MATRIX COROTATIONAL::DRDF(const MATRIX2& F, const MATRIX2& R, const MATRIX2& S, int index)
{
  MATRIX3 R3,S3;
  R3 = promote(R);
  S3 = promote(S);

  MATRIX3 eiej;
  eiej.setZero();
  const int i = index % 2;
  const int j = index / 2;
  eiej(i,j) = 1;

  MATRIX3 G = (MATRIX3::Identity() * S3.trace() - S3) * R3.transpose();
  MATRIX3 RTeiej = R3.transpose() * eiej;

  VECTOR skewVector = skew(RTeiej);
  skewVector *= 2.0;

  VECTOR omega = G.inverse() * skewVector;

  MATRIX3 cross;
  cross <<         0, -omega[2],  omega[1],
            omega[2],         0, -omega[0],
           -omega[1],  omega[0], 0;

  MATRIX3 final = cross * R3;
  return demote(final);
}

///////////////////////////////////////////////////////////////////////
// get the gradient of the rotation
///////////////////////////////////////////////////////////////////////
MATRIX COROTATIONAL::DRDF(const MATRIX2& F, int index)
{
  // get polar decomposition (this should be cached eventually)
  MATRIX2 R2,S2;
  polarDecomposition(F, R2, S2);

  MATRIX3 R,S;
  R = promote(R2);
  S = promote(S2);

  MATRIX3 eiej;
  eiej.setZero();
  const int i = index % 2;
  const int j = index / 2;
  eiej(i,j) = 1;

  MATRIX3 G = (MATRIX3::Identity() * S.trace() - S) * R.transpose();
  MATRIX3 RTeiej = R.transpose() * eiej;

  VECTOR skewVector = skew(RTeiej);
  skewVector *= 2.0;

  VECTOR omega = G.inverse() * skewVector;

  if (omega.hasNaN())
  {
    G = 2.0 * MATRIX3::Identity();
    omega = G.inverse() * skewVector;
  }

  MATRIX3 cross;
  cross <<         0, -omega[2],  omega[1],
            omega[2],         0, -omega[0],
           -omega[1],  omega[0], 0;

  /*
  cout << " RTeiej: " << endl;
  cout << RTeiej << endl;

  cout << " skew: " << endl;
  cout << skewVector << endl;

  cout << " cross: " << endl;
  cout << cross << endl;
  */

  MATRIX3 final = cross * R;
  //cout << " final: " << endl;
  //cout << final << endl;

  //return demote(cross * R);
  //return demote(final);
  //return -1.0 * demote(final);
  return demote(final);
}

///////////////////////////////////////////////////////////////////////
// promote a 2x2 matrix to 3x3, with zero off-diagonals and a
// and along the diagonal
///////////////////////////////////////////////////////////////////////
MATRIX3 COROTATIONAL::promote(const MATRIX2& A)
{
  MATRIX3 final;
  final.setZero();
  final(2,2) = 1;

  for (int y = 0; y < 2; y++)
    for (int x = 0; x < 2; x++)
      final(x,y) = A(x,y);

  return final;
}

///////////////////////////////////////////////////////////////////////
// demote a 3x3 to a 2x2, take the upper left of the matrix
///////////////////////////////////////////////////////////////////////
MATRIX2 COROTATIONAL::demote(const MATRIX3& A)
{
  MATRIX2 final;
  for (int y = 0; y < 2; y++)
    for (int x = 0; x < 2; x++)
      final(x,y) = A(x,y);

  return final;
}

///////////////////////////////////////////////////////////////////////
// nearest rotation
///////////////////////////////////////////////////////////////////////
MATRIX2 COROTATIONAL::R(const MATRIX& F)
{
  MATRIX2 R,S;
  polarDecomposition(F, R, S);
  return R;
}
  
///////////////////////////////////////////////////////////////////////
// stretching-only term
///////////////////////////////////////////////////////////////////////
MATRIX2 COROTATIONAL::S(const MATRIX& F)
{
  MATRIX2 R,scaling;
  polarDecomposition(F, R, scaling);
  return scaling;
}

///////////////////////////////////////////////////////////////////////
// derivative of rigid scaling w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX COROTATIONAL::DSDF(const MATRIX& F)
{
  MATRIX final(4,4);

  for (int i = 0; i < 4; i++)
  {
    MATRIX2 column = DSDF(F, i);
    final.col(i) = flatten(column);
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// derivative of rigid rotation w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX COROTATIONAL::DRDF(const MATRIX& F)
{
  MATRIX final(4,4);

  for (int i = 0; i < 4; i++)
  {
    MATRIX2 column = DRDF(F, i);
    final.col(i) = flatten(column);
  }

  return final;
  //return -1.0 * final;
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR COROTATIONAL::flatten(const MATRIX& A)
{
  VECTOR final(A.rows() * A.cols());

  int i = 0;
  for (int y = 0; y < A.cols(); y++)
    for (int x = 0; x < A.rows(); x++, i++)
      final[i] = A(x,y);

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the skew-symmetric portion of this matrix
///////////////////////////////////////////////////////////////////////
VECTOR COROTATIONAL::skew(const MATRIX3& A)
{
  MATRIX3 sym = (A - A.transpose()) * 0.5;

  VECTOR final(3);
  final[0] = sym(1,2);
  final[1] = sym(0,2);
  final[2] = sym(0,1);

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the gradient of the scaling
///////////////////////////////////////////////////////////////////////
MATRIX COROTATIONAL::DSDF(const MATRIX2& F, const MATRIX2& R, const MATRIX2& S, const MATRIX2& DR, int index)
{
  MATRIX2 eiej;
  eiej.setZero();
  const int i = index % 2;
  const int j = index / 2;
  eiej(i,j) = 1;

  //MATRIX2 final = R.transpose() * (eiej + DR * S);
  //MATRIX2 final = R.transpose() * eiej + DR * F;
  MATRIX2 final = R.transpose() * eiej + DR.transpose() * F;
  return final;
}

///////////////////////////////////////////////////////////////////////
// get the gradient of the scaling
///////////////////////////////////////////////////////////////////////
MATRIX COROTATIONAL::DSDF(const MATRIX2& F, int index)
#if 0
{
  MATRIX2 U,S,V;
  svd(F, U,S,V);
  MATRIX2 DU, DS, DV;
  const int i = index % 2;
  const int j = index / 2;
  svdGradient(F, i, j, DU, DS, DV);
  return DV * S * V.transpose() + V * DS * V.transpose() + V * S * DV.transpose();
}
#else
{
  MATRIX2 R,S;
  polarDecomposition(F,R,S);
  MATRIX DRDF_i = DRDF(F, index);
  
  MATRIX2 eiej;
  eiej.setZero();
  const int i = index % 2;
  const int j = index / 2;
  eiej(i,j) = 1;
  //eiej(j,i) = 1;

  //MATRIX2 final = R.transpose() * eiej + DRDF_i.transpose() * S;
  MATRIX2 final = R.transpose() * (eiej + DRDF_i * S);

  /*
  if (index < 2)
  {
    const Real a2 = final(0,1);
    final(0,1) = -a2;
    final(1,0) = -a2;
  }
  else
  {
    const Real a1 = final(1,0);
    final(0,1) = -a1;
    final(1,0) = -a1;
  }
  */

  //return final;
  return final;
  //return R.transpose() * eiej + DRDF_i.transpose() * F;
  //return R.transpose() * eiej + DRDF_i.transpose() * F;
  //return R.transpose() * eiej - R.transpose() * DRDF_i * S;
}
#endif
#if 0
{
  MATRIX2 R2,S2;
  polarDecomposition(F,R2,S2);
  //MATRIX DRDF_i = DRDF(F, index);
  //return demote(R.transpose() * (eiej - DRDF_i * S));
  
  MATRIX3 eiej;
  eiej.setZero();
  const int i = index % 2;
  const int j = index / 2;
  eiej(i,j) = 1;

  MATRIX3 R,S;
  R = promote(R2);
  S = promote(S2);

  MATRIX3 final;
  {
    // get polar decomposition (this should be cached eventually)
    MATRIX2 R2,S2;
    polarDecomposition(F, R2, S2);

    MATRIX3 R,S;
    R = promote(R2);
    S = promote(S2);

    MATRIX3 eiej;
    eiej.setZero();
    const int i = index % 2;
    const int j = index / 2;
    eiej(i,j) = 1;

    MATRIX3 G = (MATRIX3::Identity() * S.trace() - S) * R.transpose();
    MATRIX3 RTeiej = R.transpose() * eiej;

    VECTOR skewVector = skew(RTeiej);
    skewVector *= 2.0;

    VECTOR omega = G.inverse() * skewVector;

    MATRIX3 cross;
    cross <<         0, -omega[2],  omega[1],
              omega[2],         0, -omega[0],
             -omega[1],  omega[0], 0;

    final = cross * R;
  }

  //return demote(R.transpose() * (eiej - DRDF_i * S));
  return demote(R.transpose() * (eiej - final * S));
}
#endif

///////////////////////////////////////////////////////////////////////
// get the gradient of each component of the SVD
///////////////////////////////////////////////////////////////////////
void COROTATIONAL::svdGradient(const MATRIX2& F, const int i, const int j,
                                     MATRIX2& DU, MATRIX2& DS, MATRIX2& DV)
{
  MATRIX2 U,V;
  VECTOR Svec;

  svd(F, U, Svec, V);
  MATRIX2 S;
  S.setZero();
  S(0,0) = Svec[0];
  S(1,1) = Svec[1];

  DS.setZero();
  DS(0,0) = U(i, 0) * V(j, 0);
  DS(1,1) = U(i, 1) * V(j, 1);

  // solve for the off-diagonal entries of the gradients
  MATRIX2 A;
  A << Svec[0], Svec[1],
       Svec[1], Svec[0];

  MATRIX2 delta;
  delta.setZero();
  delta(i,j) = 1;

  MATRIX2 product = U.transpose() * delta * V;
  //cout << " Product matrix:" << endl;
  //cout << product << endl;
  VEC2 b;
  b[0] =  product(1,0);
  b[1] = -product(0,1);

  VEC2 x = A.inverse() * b;

  /*
  cout << " A: " << endl;
  cout << A << endl;
  cout << " b: " << endl;
  cout << b << endl;
  cout << " x: " << endl;
  cout << x << endl;
  */

  DU.setZero();
  DV.setZero();
  DU(1,0) = x[0];
  DV(1,0) = x[1];

  DU(0,1) = -x[0];
  DV(0,1) = -x[1];

  DU = U * DU;
  DV = V * DV;
  DV *= -1;
  /*
  cout << " DU: " << endl;
  cout << DU;
  cout << " DV: " << endl;
  cout << DV;
  */
}

///////////////////////////////////////////////////////////////////////
// run a unit test on the polar gradients
///////////////////////////////////////////////////////////////////////
void COROTATIONAL::testPolarGrad(const MATRIX2& F)
{
  cout << "======================================" << endl;
  cout << " Testing DRDF polar gradient " << endl;
  cout << "======================================" << endl;
  int index = 0;
  Real diff = 0;
  for (int j = 0; j < 2; j++)
    for (int i = 0; i < 2; i++, index++)
    {
      // get the DRDF
      MATRIX2 DRDF_i = DRDF(F, index);

      // get the SVD
      MATRIX2 U, V;
      VECTOR S;
      svd(F, U, S, V);

      // get the SVD gradient
      MATRIX2 DU, DS, DV;
      svdGradient(F, i, j, DU, DS, DV);

      // compare
      MATRIX2 DRDF_test = DU * V.transpose() + U * DV.transpose();
      DRDF_test *= -1;

      MATRIX diffDR = DRDF_i - DRDF_test;

#if 1
      cout << " polar version: " << endl;
      cout << DRDF_i << endl;
      cout << " SVD version: " << endl;
      cout << DRDF_test << endl;
      cout << " diff: " << endl;
      cout << diffDR << endl;
#endif

      diff += diffDR.norm();

      cout << "DRDF (" << i << ", " << j << ") diff: " << diffDR.norm() << endl;
    }

  if (diff < 1e-8)
    cout << " TEST PASSED " << endl;
  else
  {
    cout << " TEST FAILED " << endl;
    exit(0);
  }

  cout << "======================================" << endl;
  cout << " Testing DSDF polar gradient " << endl;
  cout << "======================================" << endl;
  index = 0;
  diff = 0;
  for (int j = 0; j < 2; j++)
    for (int i = 0; i < 2; i++, index++)
    {
      // get the DSDF
      MATRIX2 DSDF_i = DSDF(F, index);

      // get the SVD
      MATRIX2 U, S, V;
      svd(F, U, S, V);

      // get the SVD gradient
      MATRIX2 DU, DS, DV;
      svdGradient(F, i, j, DU, DS, DV);

      // compare
      MATRIX2 DSDF_test = DV * S * V.transpose() + V * DS * V.transpose() +
                          V * S * DV.transpose();

      MATRIX diffDS = DSDF_i - DSDF_test;
#if 0
      cout << " polar version: " << endl;
      cout << DSDF_i << endl;
      cout << " SVD version: " << endl;
      cout << DSDF_test << endl;
      cout << " diff: " << endl;
      cout << diffDS << endl;
#endif

      diff += diffDS.norm();

      cout << "DSDF (" << i << ", " << j << ") diff: " << diffDS.norm() << endl;
    }

  if (diff < 1e-8)
    cout << " TEST PASSED " << endl;
  else
  {
    cout << " TEST FAILED " << endl;
    exit(0);
  }
  //exit(0);
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
Real COROTATIONAL::psi(const MATRIX2& F)
{
  MATRIX2 R,S;
  polarDecomposition(F,R,S);

  /*
  cout << " Psi S: " << endl;
  cout << S << endl;
  cout << " Psi F: " << endl;
  cout << F << endl;
  */

  const MATRIX2 I = MATRIX2::Identity();
  //return _mu * (F - R).squaredNorm() + _lambda * 0.5 * pow((R.transpose() * F - I).trace(), 2.0);
  //return _mu * (F - R).squaredNorm() + _lambda * 0.5 * pow((R.transpose() * F - I).trace(), 2.0);
  return _mu * (F - R).squaredNorm();
}
