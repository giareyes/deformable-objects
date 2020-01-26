#include "SYMMETRIC_DIRICHLET.h"
#include <iostream>

#define USING_CONVEXIFIER 0

using namespace std;

SYMMETRIC_DIRICHLET::SYMMETRIC_DIRICHLET(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  _name = string("Symmetric Dirichlet");
  cout << " SYMMETRIC_DIRICHLET settings: " << _mu << " " << _lambda << endl;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX SYMMETRIC_DIRICHLET::PK1(const MATRIX2& F)
{
  MATRIX2 DJ;
  DJ <<  F(1,1), -F(1,0),
        -F(0,1),  F(0,0);
  Real IIS = (F.transpose() * F).trace();
  Real J = F.determinant();
  Real JSq = J * J;
  Real JCubed = JSq * J;

  return _mu * ((1.0 + 1.0 / JSq) * F - (IIS / JCubed) * DJ);
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor = 2 mu E
///////////////////////////////////////////////////////////////////////
MATRIX SYMMETRIC_DIRICHLET::PK2(const MATRIX2& F)
{
  return F.inverse() * PK1(F);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX SYMMETRIC_DIRICHLET::DPDF(const MATRIX& F)
{
  VEC4 g;
  g <<  F(1,1), -F(0,1),
        -F(1,0),  F(0,0);
  VEC4 f = flatten(F);

  const Real J = F.determinant();
  const Real J2 = J * J;
  const Real J3 = J2 * J;
  const Real J4 = J2 * J2;
  const Real IIS = (F.transpose() * F).trace();
  
  MATRIX anti(4,4);
  anti.setZero();
  anti(0,3) = 1;
  anti(1,2) = -1;
  anti(2,1) = -1;
  anti(3,0) = 1;

  MATRIX4 finalMatrix;
  finalMatrix = MATRIX4::Identity() * (1.0 + 1.0 / J2);
  finalMatrix -= (IIS / J3) * anti;
  finalMatrix += (3.0 * IIS / J4) * (g * g.transpose());
  finalMatrix -= (2.0 / J3) * (g * f.transpose() + f * g.transpose());

  Real eig2 = (J - IIS) / J3 + 1.0;
  Real eig3 = (J + IIS) / J3 + 1.0;

  // This block could be called conditionally by checking if either eig
  // is negative and then firing the SVD. However, the computation is
  // fast enough that doing it in a branch-free way here does not seem
  // too terrible.
  {
    Real projector2 = (eig2 > 0.0) ? 0.0 : eig2;
    Real projector3 = (eig3 > 0.0) ? 0.0 : eig3;
    MATRIX2 U,V;
    MATRIX2 S, middle;
    rotationVariantSVD(F, U, S, V);

    VECTOR q;
    middle << 0.0, 1.0,
             -1.0, 0.0;
    q = (1.0 / sqrt(2.0)) * flatten(U * middle * V.transpose());
    finalMatrix = finalMatrix - projector2 * (q * q.transpose());
    
    middle << 0.0, 1.0,
              1.0, 0.0;
    q = (1.0 / sqrt(2.0)) * flatten(U * middle * V.transpose());
    finalMatrix = finalMatrix - projector3 * (q * q.transpose());
  }

  // Do a sanity check, if you want.
  /*
  SelfAdjointEigenSolver<MATRIX> eigensolver(finalMatrix);
  VECTOR values = eigensolver.eigenvalues();
  for (int x = 0; x < values.size(); x++)
    if (values[x] < -1e-8)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " NEGATIVE EIGENVALUE FOUND " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Eigenvalues: \n" << values << endl;
      exit(0);
    }
    */

  return _mu * finalMatrix;
}
// A slower debugging version of the above
/*
{
  VECTOR g(4);
  g <<  F(1,1), -F(0,1),
        -F(1,0),  F(0,0);
  
  VECTOR f(4);
  f <<  F(0,0), F(1,0),
        F(0,1), F(1,1);

  const Real J = F.determinant();
  Real IIS = (F.transpose() * F).trace();
  
  MATRIX anti(4,4);
  anti.setZero();
  anti(0,3) = 1;
  anti(1,2) = -1;
  anti(2,1) = -1;
  anti(3,0) = 1;

  const Real J2 = J * J;
  const Real J3 = J2 * J;
  const Real J4 = J2 * J2;
  MATRIX finalMatrix(4,4);

  finalMatrix = MATRIX4::Identity() * (1.0 + 1.0  / J2);
  finalMatrix -= (IIS / J3) * anti;

#if 1
  finalMatrix += (3.0 * IIS / J4) * (g * g.transpose());
  finalMatrix -= (2.0 / J3) * (g * f.transpose() + f * g.transpose());
#else
  MATRIX2 U,sigma,V;
  rotationVariantSVD(F, U, sigma, V);

  // try building directly from eigensystem
  MATRIX2 DJ;
  DJ << sigma(1,1), 0,
        0, sigma(0,0);
  const Real magnitude = sqrt(sigma(0,0) * sigma(0,0) + sigma(1,1) * sigma(1,1));
  MATRIX2 Q0 = (1.0 / magnitude) * U * DJ * V.transpose();
  VECTOR q0 = flatten(Q0);
  finalMatrix += (3.0 * IIS * IIS / (J * J * J * J)) * (q0 * q0.transpose());

  // try building directly from eigensystem
  Q0 = (1.0 / sqrt(2.0)) * U * V.transpose();
  q0 = flatten(Q0);
  MATRIX2 L;
  L << 1,0,0,-1;
  MATRIX2 Q1 = (1.0 / sqrt(2.0)) * U * L * V.transpose();
  VECTOR q1 = flatten(Q1);
  const Real Jcubed = J * J * J;
  const Real summed = sigma(0,0) + sigma(1,1);
  const Real diffed = sigma(0,0) - sigma(1,1);
  finalMatrix += (2.0 / Jcubed) * (diffed * diffed * (q1 * q1.transpose()) - summed * summed * (q0 * q0.transpose()));
#endif

  SelfAdjointEigenSolver<MATRIX> eigensolverBefore(finalMatrix);
  VECTOR valuesBefore = eigensolverBefore.eigenvalues();

  Real eig2 = (J - IIS) / J3 + 1.0;
  Real eig3 = (J + IIS) / J3 + 1.0;
  Real projector2 = (eig2 > 0.0) ? 0.0 : eig2;
  Real projector3 = (eig3 > 0.0) ? 0.0 : eig3;

  MATRIX2 U,V;
  MATRIX2 S, middle;
  VECTOR q;
  rotationVariantSVD(F, U, S, V);

  middle << 0.0, 1.0,
           -1.0, 0.0;
  q = (1.0 / sqrt(2.0)) * flatten(U * middle * V.transpose());
  finalMatrix = finalMatrix - projector2 * (q * q.transpose());
  
  middle << 0.0, 1.0,
            1.0, 0.0;
  q = (1.0 / sqrt(2.0)) * flatten(U * middle * V.transpose());
  finalMatrix = finalMatrix - projector3 * (q * q.transpose());

  SelfAdjointEigenSolver<MATRIX> eigensolver(finalMatrix);
  VECTOR values = eigensolver.eigenvalues();
  for (int x = 0; x < values.size(); x++)
    if (values[x] < -1e-8)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " BAD EIGENVALUE FOUND " << endl;      
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Eigenvalues before: " << endl << valuesBefore << endl;
      cout << " Eigenvalues after: " << endl << values << endl;

      cout << " The check says (should be < 0): " << (J + J3 - IIS) << endl;
      cout << " Singular values: \n " << S << endl;
      cout << " projector 2:" << projector2 << endl;
      cout << " projector 3:" << projector3 << endl;
      exit(0);
    }

  return _mu * finalMatrix;
}
*/

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX2 SYMMETRIC_DIRICHLET::DPDFFast(const MATRIX2& F, int index)
{
  // take the polar decomposition, using SVD
  MATRIX2 R, S;
  polarDecomposition(F,R,S);

  // lots of optimizations possible here, the polar decomposition
  // and all the gradients have redundant computation
  //MATRIX2 DRDF_i = DRDF(F, index);
  MATRIX2 DRDF_i = DRDF(F, R, S, index);
  MATRIX2 DSDF_i = DSDF(F, R, S, DRDF_i, index);

  // the necessity for this is a little mysterious. The finite difference
  // verification for DRDF does not bear out the negative, but it is needed
  // for the DPDF and K verifications.
  DRDF_i *= -1;

  const MATRIX2 I = MATRIX2::Identity();
  MATRIX2 final;
  final = DRDF_i * (2.0 * _mu * (S - I));
  final += R * (2.0 * _mu * DSDF_i);
  
  return final;
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR SYMMETRIC_DIRICHLET::flatten(const MATRIX& A)
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
Real SYMMETRIC_DIRICHLET::psi(const MATRIX2& F)
{
  /*
  Real IIIS = F.determinant();

  MATRIX2 R,S;
  polarDecomposition(F,R,S);
  Real IS = S.trace();

  Real IIS = (F.transpose() * F).trace();

#if USING_CONVEXIFIER
  if (IS < 2.0)
    return 0.5 * _mu * (F + R).squaredNorm();
#endif

  return 0.5 * _mu * (F - R).squaredNorm();
  */
  Real IIS = (F.transpose() * F).trace();
  Real J = F.determinant();

  return _mu * 0.5 * (IIS + IIS / (J * J));
}

///////////////////////////////////////////////////////////////////////
// take the svd of F, so we can form the polar decomposition
///////////////////////////////////////////////////////////////////////
void SYMMETRIC_DIRICHLET::svd(const MATRIX2& F, MATRIX2& U, VECTOR& S, MATRIX2& V) const
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
bool SYMMETRIC_DIRICHLET::symSchur2(const MATRIX2& A, Real& c, Real& s) const
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
    c = 1;
    s = 0;
    return false;
  }
}

///////////////////////////////////////////////////////////////////////
// take the svd of F, so we can form the polar decomposition
///////////////////////////////////////////////////////////////////////
void SYMMETRIC_DIRICHLET::mySVD(const MATRIX2& F, MATRIX2& U, VECTOR& S, MATRIX2& V) const
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
void SYMMETRIC_DIRICHLET::rotationVariantSVD(const MATRIX2& F, MATRIX2& U, MATRIX2& sigma, MATRIX2& V) const
{
  VECTOR s;
  svd(F, U,s,V);

  // pipe the singular values to a matrix
  sigma << s[0], 0,
           0,    s[1];

  MATRIX2 Uoriginal = U;
  MATRIX2 sigmaOriginal = sigma;

  // reflection matrix
  MATRIX2 L;
  L.setZero();
  L(0,0) = 1;
  L(1,1) = (U * V.transpose()).determinant();

  const Real detU = U.determinant();
  const Real detV = V.determinant();

  // which to pull the reflection out of?
  if (detU < 0 && detV > 0)
    U = U * L;
  else if (detU > 0 && detV < 0)
    V = V * L;
  sigma = sigma * L;
}

///////////////////////////////////////////////////////////////////////
// take the polar decomposition, using SVD
///////////////////////////////////////////////////////////////////////
void SYMMETRIC_DIRICHLET::polarDecomposition(const MATRIX2& F, MATRIX2& R, MATRIX2& S) const
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
  L(0,0) = 1;
  L(1,1) = (U * V.transpose()).determinant();

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
MATRIX SYMMETRIC_DIRICHLET::DSDF(const MATRIX2& F, const MATRIX2& R, const MATRIX2& S, const MATRIX2& DR, int index)
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
MATRIX3 SYMMETRIC_DIRICHLET::promote(const MATRIX2& A)
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
MATRIX2 SYMMETRIC_DIRICHLET::demote(const MATRIX3& A)
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
VECTOR SYMMETRIC_DIRICHLET::skew(const MATRIX3& A)
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
MATRIX SYMMETRIC_DIRICHLET::DRDF(const MATRIX2& F, const MATRIX2& R, const MATRIX2& S, int index)
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
Real SYMMETRIC_DIRICHLET::DpsiDF(const MATRIX2& F, const int index)
{
  // take the polar decomposition, using SVD
  MATRIX2 R, S;
  polarDecomposition(F,R,S);
  
  // lots of optimizations possible here, the polar decomposition
  // and all the gradients have redundant computation
  MATRIX2 DRDF_i = DRDF(F, R, S, index);
  MATRIX2 DSDF_i = DSDF(F, R, S, DRDF_i, index);

  MATRIX2 left = DSDF_i - MATRIX2::Identity();
  MATRIX2 right = S - MATRIX2::Identity();

  Real final = 0;
  for (int x = 0; x < 2; x++)
    for (int y = 0; y < 2; y++)
      final += left(x,y) * right(x,y);

  return 2 * final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
bool SYMMETRIC_DIRICHLET::testSVD()
{
  MATRIX2 F;
  F << 0.710625, 0.00236432,
      0.0687256,    1.00035;

  SYMMETRIC_DIRICHLET length;

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

///////////////////////////////////////////////////////////////////////
// solve for the roots of a quadratic
// should probably put some guards in here later on
///////////////////////////////////////////////////////////////////////
VEC2 quadraticRoots(const VEC3& coeffs)
{
  const Real& a = coeffs[0];
  const Real& b = coeffs[1];
  const Real& c = coeffs[2];

  const Real inv2a = 1.0 / (2.0 * a);
  const Real insideRadical = b * b - 4.0 * a * c;
  const Real radical = sqrt(insideRadical) * inv2a;
  const Real front = -b * inv2a;
  VEC2 roots;
  roots[0] = front + radical;
  roots[1] = front - radical;

  return roots; 
}

///////////////////////////////////////////////////////////////////////
// get the eigenvalues of DPDF, if we have expressions for them
///////////////////////////////////////////////////////////////////////
VEC4 SYMMETRIC_DIRICHLET::eigenvalues(const MATRIX2& F)
{
  Real IC = (F.transpose() * F).trace();
  Real J = F.determinant();
  Real J2 = J * J;
  Real J3 = J2 * J;

  //Real J4 = J2 * J2;
  //Real J6 = J4 * J2;

  // note that since singular values are all taken to an
  // even power, the rotation variant version is not needed
  MATRIX2 U,V;
  VECTOR S;
  //rotationVariantSVD(F, U, sigma, V);
  svd(F, U, S, V);
  Real s0 = S[0];
  Real s1 = S[1];

  /*
  // solve the quadratic the direct way
  VEC3 coeffs;
  coeffs[0] = 1.0;
  coeffs[1] = -3.0 * IC * IC / J4 + 8.0 / J2;
  coeffs[2] = -3.0 * IC * IC / J6 + 16.0 / J4;
  VEC2 quadratic = quadraticRoots(coeffs);

  cout << " roots: " << quadratic << endl;

  VEC2 direct;
  direct[0] = 3.0 / pow(s0, 4.0) - 1.0 / J2;
  direct[1] = 3.0 / pow(s1, 4.0) - 1.0 / J2;
  */

  const Real s0_4 = s0 * s0 * s0 * s0;
  const Real s1_4 = s1 * s1 * s1 * s1;
  VEC4 values;
  values[0] = _mu * (3.0 / s0_4 + 1.0);
  values[1] = _mu * (3.0 / s1_4 + 1.0);
  values[2] = _mu * ((J + IC) / J3 + 1.0);
  values[3] = _mu * ((J - IC) / J3 + 1.0);

  std::sort(values.data(),values.data()+values.size());

  return values;
}

static void sortEigensystem(MATRIX4& vectors, VEC4& values)
{
  vector<pair<Real,int> > sortedPairs;
  for (int i = 0; i < 4; i++)
  {
    pair<Real,int> P = make_pair(values[i], i);
    sortedPairs.push_back(P);
  }
  sort(sortedPairs.begin(), sortedPairs.end());

  cout << "Sorted indices: " << endl;
  for (int i = 0 ; i < 4 ; i++) 
    cout << sortedPairs[i].second << endl;

  VEC4 sortedValues;
  MATRIX4 sortedVectors;

  for (int i = 0; i < 4; i++)
  {
    int index = sortedPairs[i].second;
    sortedValues[i] = values[index];
    sortedVectors.col(i) = vectors.col(index);
  }

  values = sortedValues;
  vectors = sortedVectors;
}

///////////////////////////////////////////////////////////////////////
// get the eigensysem of DPDF, if we have expressions for them
///////////////////////////////////////////////////////////////////////
void SYMMETRIC_DIRICHLET::eigensystem(const MATRIX2& F, MATRIX4& vectors, VEC4& values)
{
  Real IC = (F.transpose() * F).trace();
  Real J = F.determinant();
  Real J2 = J * J;
  Real J3 = J2 * J;

  // note that since singular values are all taken to an
  // even power, the rotation variant version is not needed,
  //
  // however we use it here because we still need to pull
  // the reflection off of the U or V. The convention doesn't
  // matter though; it doesn't care which singular value the
  // minus is then loaded into.
  MATRIX2 U,V;
  MATRIX2 S;
  rotationVariantSVD(F, U, S, V);
  Real s0 = S(0,0);
  Real s1 = S(1,1);

  const Real s0_4 = s0 * s0 * s0 * s0;
  const Real s1_4 = s1 * s1 * s1 * s1;
  values[0] = _mu * (3.0 / s0_4 + 1.0);
  values[1] = _mu * (3.0 / s1_4 + 1.0);
  values[2] = _mu * ((J + IC) / J3 + 1.0);
  values[3] = _mu * ((J - IC) / J3 + 1.0);

  // some temp variables for constructing the eigensystems
  MATRIX2 middle;
  MATRIX2 Q;
  VEC4 q;

  // first eigenvector
  middle << 1.0, 0.0,
            0.0, 0.0;
  Q = U * middle * V.transpose();
  q = flatten(Q);
  vectors.col(0) = q;

  // second eigenvector
  middle << 0.0, 0.0,
            0.0, 1.0;
  Q = U * middle * V.transpose();
  q = flatten(Q);
  vectors.col(1) = q;

  // third eigenvector
  middle << 0.0, 1.0,
            1.0, 0.0;
  Q = U * middle * V.transpose();
  q = flatten(Q);
  vectors.col(2) = (1.0 / sqrt(2.0)) * q;

  // fourth eigenvector
  middle << 0.0,  1.0,
            -1.0, 0.0;
  Q = U * middle * V.transpose();
  q = flatten(Q);
  vectors.col(3) = (1.0 / sqrt(2.0)) * q;

  sortEigensystem(vectors, values);
}
