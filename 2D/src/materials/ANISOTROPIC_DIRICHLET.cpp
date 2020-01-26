#include "ANISOTROPIC_DIRICHLET.h"
#include <iostream>

using namespace std;

#define USING_SIGNED_VERSION 1

ANISOTROPIC_DIRICHLET::ANISOTROPIC_DIRICHLET(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  _name = string("ANISOTROPIC_DIRICHLET");
  //cout << " ANISOTROPIC_DIRICHLET settings are mu = " << mu << endl;

  _a = VEC2(1.0, 1.0);
  _a *= 1.0 / _a.norm();

  //cout << " Anisotropy direction: " << _a << endl;

  _A = _a * _a.transpose();
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX ANISOTROPIC_DIRICHLET::PK1(const MATRIX2& F)
{
  return _mu * 2.0 * F * _A;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void ANISOTROPIC_DIRICHLET::setAnisotropyDirection(const VEC2& a)
{
  _a = a;
  _a *= 1.0 / _a.norm();
  _A = _a * _a.transpose();
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor = 2 mu E
///////////////////////////////////////////////////////////////////////
MATRIX ANISOTROPIC_DIRICHLET::PK2(const MATRIX2& F)
{
  return F.inverse() * PK1(F);
}

static VECTOR vec(const MATRIX2 A)
{
  VECTOR v(4);
  v[0] = A(0,0);
  v[1] = A(1,0);
  v[2] = A(0,1);
  v[3] = A(1,1);

  return v;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX ANISOTROPIC_DIRICHLET::DPDF(const MATRIX& F)
{
  MATRIX H(4,4);

  MATRIX blockA(4,4);
  blockA.setZero();
  blockA(0,0) = _A(0,0);
  blockA(1,1) = _A(0,0);

  blockA(2,0) = blockA(0,2) = _A(1,0);
  blockA(3,1) = blockA(1,3) = _A(1,0);

  blockA(2,2) = _A(1,1);
  blockA(3,3) = _A(1,1);
  H = 2.0 * _mu * blockA;
  
  return H;
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR ANISOTROPIC_DIRICHLET::flatten(const MATRIX& A)
{
  VECTOR final(A.rows() * A.cols());

  int i = 0;
  for (int y = 0; y < A.cols(); y++)
    for (int x = 0; x < A.rows(); x++, i++)
      final[i] = A(x,y);

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
Real ANISOTROPIC_DIRICHLET::psi(const MATRIX2& F)
{
  Real IV = _a.transpose() * (F.transpose() * F) * _a;
  return _mu * IV;
}

///////////////////////////////////////////////////////////////////////
// take the svd of F, so we can form the polar decomposition
///////////////////////////////////////////////////////////////////////
void ANISOTROPIC_DIRICHLET::svd(const MATRIX2& F, MATRIX2& U, VECTOR& S, MATRIX2& V) const
{
  JacobiSVD<MATRIX> svd(F, ComputeFullU | ComputeFullV);
  U = svd.matrixU();
  V = svd.matrixV();
  S = svd.singularValues();
}

///////////////////////////////////////////////////////////////////////
// generate coefficients for Jacobi (Givens) rotation
// Algorithm 8.4.1 from Golub and Van Loan 1996
///////////////////////////////////////////////////////////////////////
bool ANISOTROPIC_DIRICHLET::symSchur2(const MATRIX2& A, Real& c, Real& s) const
{
  //const Real eps = (Real)2 * std::numeric_limits<Real>::epsilon();
  //const Real eps = (Real)std::numeric_limits<Real>::epsilon();
  const Real eps = (Real)std::numeric_limits<float>::epsilon();

  const int p = 1;
  const int q = 0;

  cout << " A(p,q): " << fabs(A(p,q)) << endl;
  cout << " eps * sqrt(A(q,q) * A(p,p)): " << fabs(eps * sqrt(A(q,q) * A(p,p))) << endl;

  if (fabs(A(p,q)) > fabs(eps * sqrt(A(q,q) * A(p,p))))
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    Real tau = (A(q,q) - A(p,p)) / (2.0 * A(p,q));

    /*
    Real sign = (tau > 0) ? 1 : -1;
    Real t = sign  / (fabs(tau) + sqrt(1.0 + tau * tau));
    */

    Real t;
    if (tau >= 0.0)
      t = 1.0 / (tau + sqrt(1.0 + tau * tau));
    else
      t = 1.0 / (tau - sqrt(1.0 + tau * tau));

    c = 1.0 / sqrt(1.0 + t * t);
    s = t * c;
    return true;
  }
  else
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    c = 1;
    s = 0;
    return false;
  }
}

///////////////////////////////////////////////////////////////////////
// take the svd of F, so we can form the polar decomposition
///////////////////////////////////////////////////////////////////////
void ANISOTROPIC_DIRICHLET::mySVD(const MATRIX2& F, MATRIX2& U, VECTOR& S, MATRIX2& V) const
{
  // form A^T * A
  MATRIX2 A = F.transpose() * F;

  // init V before the rotation
  V = MATRIX2::Identity();

  int sweeps = 0;

  const Real considerAsZero = (Real)2 * std::numeric_limits<Real>::denorm_min();
  //const Real precision = (Real)2 * std::numeric_limits<Real>::epsilon();
  const Real precision = std::numeric_limits<float>::epsilon();

  Real maxDiag = (fabs(A(0,0)) > fabs(A(1,1))) ? precision * fabs(A(0,0)) : precision * fabs(A(1,1));
  Real threshold = (considerAsZero > maxDiag) ? considerAsZero : maxDiag;

  //Real old = fabs(2.0 * A(1,0));
  //while (fabs(A(1,0)) < old)
  bool success = true;
  //while (fabs(A(1,0)) > fabs(precision * sqrt(A(0,0) * A(1,1))))
  while (success)
  {
    cout << "========================================" << endl; 
    cout << " sweep: " << sweeps << endl;
    cout << " F: " << endl << F << endl;
    cout << " A: " << endl << A << endl;
    Real before = fabs(A(1,0));
    cout << " A(1,0) before: " << fabs(A(1,0)) << endl;

    // kill off the one diagonal
    Real c,s;
    success = symSchur2(A, c, s);

    if (success)
    {
      // build the Givens rotation
      MATRIX2 J;
      J <<  c, -s,
            s,  c;

      // annihilate the off-diagonal term
      V = V * J;
      A = J.transpose() * A * J;
      //A = J * A * J.transpose();
      sweeps++;
     
      cout << " J: " << endl << J << endl;
      cout << " V: " << endl << V << endl;
      cout << " c,s: " << c << ", " << s << endl;
      cout << " A after: " << endl << A << endl;
      Real after = fabs(A(1,0));
      cout << " A(1,0) after: " << fabs(A(1,0)) << endl;
    }

    /*
    if (after > before && after > 1e-1)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " JACOBI ROTATION FAILED " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      exit(0);
    }
    */

    maxDiag = (fabs(A(0,0)) > fabs(A(1,1))) ? precision * fabs(A(0,0)) : fabs(A(1,1));
    threshold = (considerAsZero > maxDiag) ? considerAsZero : maxDiag;

    if (sweeps > 5)
      exit(0);
  }
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " Sweeps: " << sweeps << endl;

  // form the sigma
  S.resize(2);
  S[0] = sqrt(A(0,0));
  S[1] = sqrt(A(1,1));

  // see if the columns need a swap
  if (S[0] < S[1])
  {
    Real temp = S[0];
    S[0] = S[1];
    S[1] = temp;
    VEC2 column = V.col(0);
    V.col(0) = V.col(1);
    V.col(1) = column;
  }

  U = F * V;
  U.col(0) *= 1.0 / S[0];
  U.col(1) *= 1.0 / S[1];

  if (fabs(S[1]) < precision)
  {
    VEC2 col = U.col(0);
    VEC2 reflect;
    //reflect << col[1], -col[0];
    reflect << -col[1], col[0];
    U.col(1) = reflect;
  }

  // DEBUG: sanity check
  MATRIX2 Sigma;
  Sigma << S[0], 0,
              0, S[1];
  MATRIX2 test = U * Sigma * V.transpose();

  if ((F - test).norm() > 2 * precision)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " BAD SVD " << endl;
    cout << " F: " << endl << F << endl;
    cout << " reconstructed: " << endl << test << endl;
    cout << " Diff: " << F - test << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////
// take the polar decomposition, using SVD
///////////////////////////////////////////////////////////////////////
void ANISOTROPIC_DIRICHLET::polarDecomposition(const MATRIX2& F, MATRIX2& R, MATRIX2& S) const
{
  MATRIX2 U;
  VECTOR s;
  MATRIX2 V;
  svd(F, U,s,V);

  // pipe the singular values to a matrix
  MATRIX2 sigma;
  sigma << s[0], 0,
           0,    s[1];

  MATRIX2 Uoriginal = U;
  MATRIX2 sigmaOriginal = sigma;

  // reflection matrix
  MATRIX2 L;
  L.setZero();

  Real product = (U * V.transpose()).determinant();

  L(0,0) = 1;
  L(1,1) = product;
  /*
  if (product < 0.0 && sigma.trace() > 2.0)
  {
    L(0,0) = 1;
    L(1,1) = product;
  }
  else
  {
    L(0,0) = product;
    L(1,1) = 1;
  }
  */

  const Real detU = U.determinant();
  const Real detV = V.determinant();

  // which to pull the reflection out of?
  if (detU < 0 && detV > 0)
    U = U * L;
  else if (detU > 0 && detV < 0)
    V = V * L;
  sigma = sigma * L;

  R = U * V.transpose();
  S = V * sigma * V.transpose();
}

///////////////////////////////////////////////////////////////////////
// get the gradient of the scaling
///////////////////////////////////////////////////////////////////////
MATRIX ANISOTROPIC_DIRICHLET::DSDF(const MATRIX2& F, const MATRIX2& R, const MATRIX2& S, const MATRIX2& DR, int index)
{
  MATRIX2 eiej;
  eiej.setZero();
  const int i = index % 2;
  const int j = index / 2;
  eiej(i,j) = 1;

  MATRIX2 final = R.transpose() * (eiej + DR * S);
  return final;
}

///////////////////////////////////////////////////////////////////////
// promote a 2x2 matrix to 3x3, with zero off-diagonals and a
// and along the diagonal
///////////////////////////////////////////////////////////////////////
MATRIX3 ANISOTROPIC_DIRICHLET::promote(const MATRIX2& A)
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
MATRIX2 ANISOTROPIC_DIRICHLET::demote(const MATRIX3& A)
{
  MATRIX2 final;
  for (int y = 0; y < 2; y++)
    for (int x = 0; x < 2; x++)
      final(x,y) = A(x,y);

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the skew-symmetric portion of this matrix
///////////////////////////////////////////////////////////////////////
VECTOR ANISOTROPIC_DIRICHLET::skew(const MATRIX3& A)
{
  MATRIX3 sym = (A - A.transpose()) * 0.5;

  VECTOR final(3);
  final[0] = sym(1,2);
  final[1] = sym(0,2);
  final[2] = sym(0,1);

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the gradient of the rotation
///////////////////////////////////////////////////////////////////////
MATRIX ANISOTROPIC_DIRICHLET::DRDF(const MATRIX2& F, const MATRIX2& R, const MATRIX2& S, int index)
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

  if (omega.hasNaN())
  {
    G = 2.0 * MATRIX3::Identity();
    omega = G.inverse() * skewVector;
  }

  /*
  cout << " S3: " << endl << S3 << endl;
  cout << " S3.trace(): " << S3.trace() << endl;
  cout << " G: " << endl << G << endl;
  cout << " inverse: " << G.inverse() << endl;
  */

  MATRIX3 cross;
  cross <<         0, -omega[2],  omega[1],
            omega[2],         0, -omega[0],
           -omega[1],  omega[0], 0;

  MATRIX3 final = cross * R3;
  return demote(final);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
bool ANISOTROPIC_DIRICHLET::testSVD()
{
  MATRIX2 F;
  F << 0.710625, 0.00236432,
      0.0687256,    1.00035;

  ANISOTROPIC_DIRICHLET length;

  MATRIX2 U,V;
  VECTOR S;
  length.mySVD(F, U, S, V);
  cout << " Mine ============================ " << endl;
  cout << "U: " << endl << U << endl;
  cout << "S: " << endl << S << endl;
  cout << "V: " << endl << V << endl;
  length.svd(F, U, S, V);
  cout << " Eigen ============================ " << endl;
  cout << "U: " << endl << U << endl;
  cout << "S: " << endl << S << endl;
  cout << "V: " << endl << V << endl;

  F << 0, 0,
       0, 1;

  length.mySVD(F, U, S, V);
  cout << " Mine ============================ " << endl;
  cout << "U: " << endl << U << endl;
  cout << "S: " << endl << S << endl;
  cout << "V: " << endl << V << endl;
  length.svd(F, U, S, V);
  cout << " Eigen ============================ " << endl;
  cout << "U: " << endl << U << endl;
  cout << "S: " << endl << S << endl;
  cout << "V: " << endl << V << endl;

  return false;
}
