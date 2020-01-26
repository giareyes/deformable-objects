#include "SYMMETRIC_ARAP.h"
#include <iostream>

using namespace std;

SYMMETRIC_ARAP::SYMMETRIC_ARAP(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  _name = string("Symmetric ARAP");
  cout << " SYMMETRIC_ARAP settings: " << _mu << " " << _lambda << endl;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX SYMMETRIC_ARAP::PK1(const MATRIX2& F)
{
  MATRIX2 R,S;
  polarDecomposition(F,R,S);
  
  MATRIX2 U,sigma,V;
  rotationVariantSVD(F, U, sigma, V);

  MATRIX2 DJ;
  DJ <<  F(1,1), -F(1,0),
        -F(0,1),  F(0,0);
  Real IS   = S.trace();
  Real IIS  = (F.transpose() * F).trace();
  Real IIIS = F.determinant();
  Real IIIS2 = IIIS * IIIS;

  MATRIX2 twist;
  twist << 0, 1,
           -1, 0;
  MATRIX2 Q = U * twist * V.transpose();

  MATRIX2 P = (1.0 + 1.0 / IIIS2) * F - 
              (1.0 + 1.0 / IIIS) * R +
              (1.0 / IIIS2) * (IS - IIS / IIIS) * DJ -
              (1.0 / IIIS) * (2.0 / IS) * (Q * DJ.transpose()).trace() * Q;

  return _mu * P;
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor = 2 mu E
///////////////////////////////////////////////////////////////////////
MATRIX SYMMETRIC_ARAP::PK2(const MATRIX2& F)
{
  return F.inverse() * PK1(F);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX SYMMETRIC_ARAP::DPDF(const MATRIX& F)
{
  MATRIX4 vectors;
  VEC4 values;
  eigensystem(F, vectors, values);

  // do the clamping here
  for (int x = 0; x < 4; x++)
    values[x] = (values[x] < 0.0) ? 0.0 : values[x];

  MATRIX4 finalMatrix;
  finalMatrix = vectors * values.asDiagonal() * vectors.transpose();
  finalMatrix = _mu * finalMatrix;

  // Do a sanity check, if you want.
  SelfAdjointEigenSolver<MATRIX> eigensolver(finalMatrix);
  VECTOR numerical = eigensolver.eigenvalues();
  for (int x = 0; x < values.size(); x++)
    if (numerical[x] < -1e-7)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " NEGATIVE EIGENVALUE FOUND " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Eigenvalues: \n" << numerical << endl;

      cout << " The direct ones were: \n" << values << endl;
      exit(0);
    }

  return finalMatrix;
}
// this is the direct construction version, and its correctness has been
// verified. The eigenvector outer product version above is probably the 
// preferred version for projected Newton
/*
{
  VEC4 g;
  g <<  F(1,1), -F(0,1),
        -F(1,0),  F(0,0);
  VEC4 f = flatten(F);

  MATRIX2 U,V, sigma;
  rotationVariantSVD(F, U, sigma, V);
  VECTOR q;
  MATRIX2 middle;
  middle << 0.0, 1.0,
           -1.0, 0.0;
  q = (1.0 / sqrt(2.0)) * flatten(U * middle * V.transpose());
  VECTOR r(4);
  r = flatten(U * V.transpose());

  const Real J = F.determinant();
  const Real IS = sigma.trace();
  const Real IIS = (F.transpose() * F).trace();
  const Real J2 = J * J;
  const Real J3 = J * J * J;
  
  MATRIX anti(4,4);
  anti.setZero();
  anti(0,3) = 1;
  anti(1,2) = -1;
  anti(2,1) = -1;
  anti(3,0) = 1;

  MATRIX4 finalMatrix;

  finalMatrix = MATRIX4::Identity() * (1.0 + 1.0 / J2);
  finalMatrix -= (2.0 / J3) * (f * g.transpose() + g * f.transpose());
  finalMatrix += (1.0 / J2) * (r * g.transpose() + g * r.transpose());
  finalMatrix += (1.0 / J2) * (IS - IIS / J) * anti;
  finalMatrix += (1.0 / J3) * (3.0 * IIS / J - 2.0 * IS) * (g * g.transpose());
  finalMatrix -= (1.0 + 1.0 / J) * (2.0 / IS) * (q * q.transpose());

  // Do a sanity check, if you want.
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

  return _mu * finalMatrix;
}
*/

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR SYMMETRIC_ARAP::flatten(const MATRIX& A)
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
Real SYMMETRIC_ARAP::psi(const MATRIX2& F)
{
  MATRIX2 R,S;
  polarDecomposition(F,R,S);

  //return 0.5 * _mu * ((F + R).squaredNorm() + (F.inverse() - R.transpose()).squaredNorm());

  // using this form, because otherwise some R' * R products cause
  // irritating constants to appear
  Real IS = S.trace();
  Real IIS = (F.transpose() * F).trace();
  Real J = F.determinant();

  return 0.5 * _mu * (IIS - 2.0 * IS + IIS / (J * J) - 2 * IS / J);
}

///////////////////////////////////////////////////////////////////////
// take the svd of F, so we can form the polar decomposition
///////////////////////////////////////////////////////////////////////
void SYMMETRIC_ARAP::svd(const MATRIX2& F, MATRIX2& U, VECTOR& S, MATRIX2& V) const
{
  JacobiSVD<MATRIX> svd(F, ComputeFullU | ComputeFullV);
  U = svd.matrixU();
  V = svd.matrixV();
  S = svd.singularValues();
}

///////////////////////////////////////////////////////////////////////
// take the polar decomposition, using SVD
///////////////////////////////////////////////////////////////////////
void SYMMETRIC_ARAP::rotationVariantSVD(const MATRIX2& F, MATRIX2& U, MATRIX2& sigma, MATRIX2& V) const
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
void SYMMETRIC_ARAP::polarDecomposition(const MATRIX2& F, MATRIX2& R, MATRIX2& S) const
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
// get the eigenvalues of DPDF, if we have expressions for them
///////////////////////////////////////////////////////////////////////
VEC4 SYMMETRIC_ARAP::eigenvalues(const MATRIX2& F)
{
  MATRIX2 U,sigma,V;
  rotationVariantSVD(F, U, sigma, V);
  const Real s0 = sigma(0,0);
  const Real s1 = sigma(1,1);
  const Real IS = sigma.trace();
  const Real IIS = (F.transpose() * F).trace();
  const Real J = F.determinant();
  const Real J2 = J * J;
  const Real J3 = J2 * J;

  const Real s0_3 = s0 * s0 * s0;
  const Real s1_3 = s1 * s1 * s1;
  const Real s0_4 = s0_3 * s0;
  const Real s1_4 = s1_3 * s1;
  VEC4 values;
  values[0] = _mu * (1.0 - 2.0 / s0_3 + 3.0 / s0_4);
  values[1] = _mu * (1.0 - 2.0 / s1_3 + 3.0 / s1_4);
  values[2] = _mu * (1.0 + (1.0 / J2) * (IIS / J - IS + 1));
  values[3] = _mu * (1.0 + (1.0 / IS) * ((IIS + IS) / J2 - (IIS * IS) / J3 - 2.0));

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

  /*
  cout << "Sorted indices: " << endl;
  for (int i = 0 ; i < 4 ; i++) 
    cout << sortedPairs[i].second << endl;
    */

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
void SYMMETRIC_ARAP::eigensystem(const MATRIX2& F, MATRIX4& vectors, VEC4& values)
{
  MATRIX2 U,sigma,V;
  rotationVariantSVD(F, U, sigma, V);
  const Real s0 = sigma(0,0);
  const Real s1 = sigma(1,1);
  const Real IS = sigma.trace();
  const Real IIS = (F.transpose() * F).trace();
  const Real J = F.determinant();
  const Real J2 = J * J;
  const Real J3 = J2 * J;

  const Real s0_3 = s0 * s0 * s0;
  const Real s1_3 = s1 * s1 * s1;
  const Real s0_4 = s0_3 * s0;
  const Real s1_4 = s1_3 * s1;
  values[0] = _mu * (1.0 - 2.0 / s0_3 + 3.0 / s0_4);
  values[1] = _mu * (1.0 - 2.0 / s1_3 + 3.0 / s1_4);
  values[2] = _mu * (1.0 + (1.0 / J2) * (IIS / J - IS + 1));
  values[3] = _mu * (1.0 + (1.0 / IS) * ((IIS + IS) / J2 - (IIS * IS) / J3 - 2.0));

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

  //sortEigensystem(vectors, values);
}
