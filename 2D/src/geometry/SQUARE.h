#ifndef SQUARE_H
#define SQUARE_H

#include "SETTINGS.h"
#include "MATERIAL.h"
#include "COROTATIONAL.h"
#include <vector>

class SQUARE 
{
public:
  SQUARE(MATERIAL* material, const std::vector<VEC2*>& vertices);

  void draw();
  void drawForces();
  void drawInversionDirection();

  // get the deformation gradient at Xi, the normalized position
  // inside a [-1,-1], [1,1] square
  MATRIX2 computeF(const VEC2& Xi) const;

  // get the average F at all the Gauss points
  MATRIX2 gaussF();

  void scale(const VEC2& scalar);
  void translate(const VEC2& translation);
  void rotate(const Real angle);

  VEC2 vertexAverage();

  void computeForces();
  VECTOR computeForceVector();
  VECTOR computeForceVectorFast();
  VECTOR computeInvertibleForceVector();
  MATRIX computeForceJacobian();
  MATRIX computeForceJacobianFast();
  MATRIX computeInvertibleForceJacobian();
  MATRIX computeDPDu();

  MATRIX secondPiolaKirchhoff();

  const MATRIX forces() { return _forces; };
  VECTOR forceVector() const;
  Real& conditionNumber() { return _conditionNumber; };

  VECTOR getDisplacement() const;
  void setDisplacement(VECTOR& u);

  // get pointer to the vertex
  VEC2* vertex(int i) { return _vertices[i]; };

  // get the area of the square
  Real area() const;

  // get the strain energy inside the square
  Real psi() const;
  
  // run some SVD-related tests
  void svdTest();

  // update the inversion direction information
  void updateInversionDirections();

  void printReflectionInformation();

  // run all the unit tests
  static void runUnitTests(MATERIAL* material);  
  static void runCorotationalUnitTests(COROTATIONAL* corotational);  
  static void runHypertanUnitTests(MATERIAL* material);  
  static void runEigenvalueTest(MATERIAL* material);  

private:
  // pack the vertices into a matrix
  MATRIX vertexMatrix() const;

  // get the shape function derivative matrix H for computing F
  MATRIX computeH(const VEC2& Xi) const;

  // get the forces at a Gauss point
  MATRIX computeForces(const VEC2& gaussPoint);
  VECTOR computeForceVector(const VEC2& gaussPoint);
  VECTOR computeForceVectorFast(const int gaussPoint);
  VECTOR computeInvertibleForceVector(const int gaussPoint);
  
  // get the forces Jacobian at a Gauss point
  MATRIX computeForceJacobian(const VEC2& gaussPoint);
  MATRIX computeForceJacobianFast(const int gaussPoint);
  MATRIX computeInvertibleForceJacobian(const int gaussPoint);

  MATRIX computeDPDu(const VEC2& gaussPoint, const int column);
  MATRIX computeDSDu(const VEC2& gaussPoint, const int column);

  // sample Piola-Kirchoff for unit testing
  MATRIX firstPiolaKirchhoff(const VEC2& gaussPoint);
  MATRIX secondPiolaKirchhoff(const VEC2& gaussPoint);

  // compute the derivative of deformation gradient F
  // with respect to displacement u
  MATRIX computeDFDu(const VEC2& Xi);
  MATRIX computeDFDu(const VEC2& Xi, int index);
  MATRIX computeDFDu(const int gaussPoint, int index);

  // compute the derivative of E = 1/2 (F^T * F - I)
  // with respect to F
  MATRIX computeDEDF(const MATRIX2& F);

  // compute the derivative of S, the second Piola-Kirchhoff
  // with respect to E = 1/2 (F^T * F - I)
  MATRIX DSDE();
  
  // flatten a matrix into a vector, stacking each of the columns
  // on top of each other
  //
  // this should be moved to MATERIAL
  static VECTOR flatten(const MATRIX& A);
  
  // unflatten a vector into a matrix, stacking each of the columns
  // on top of each other
  static MATRIX unflatten(const VECTOR& v, const int rows = 2);

  // take the entries of the matrix and put them along the diagonal
  MATRIX diag(const MATRIX& A);

  // take the gradient of FS (i.e. PK1) w.r.t F assuming that S is frozen
  // \frac{\partial F S}{\partial F}
  MATRIX DFSDF(const MATRIX& S);
  
  // repeat the given matrix "repeats" times along the block diagonal
  MATRIX blockDiag(const MATRIX& A, const int repeats);

  // color ramp function (used to show the condition number)
  VEC3 ramp(float value);

  ////////////////////////////////////////////////////////////////////////////
  // Inversion handling
  ////////////////////////////////////////////////////////////////////////////
  bool _oldClamped[4];
  VEC2 _oldDirection[4];
  bool _newClamped[4];
  VEC2 _newDirection[4];
  VEC2 _drawDirection[4];

  ////////////////////////////////////////////////////////////////////////////
  // UNIT TESTS
  ////////////////////////////////////////////////////////////////////////////
  void squash();
  void testK();
  void testDPDu();
  void testDPDF();
  void testDSDu();
  void testDFDu();
  void testDEDF();
  void testDEDu();
  void testDpsiDF();
 
  // see if the forces disappear at I
  void testRestStability();

  void testDRDF();  // only meaningful for COROTATIONAL material
  void testDSDF();  // only meaningful for COROTATIONAL material
  void testSVD();   // only meaningful for COROTATIONAL material
  void testMu();   // only meaningful for HYPER_TAN material
  
  // see if eigenvalues match up to the call in "eigenvalues"
  void testEigenvalues();
  
  // see if eigenvalues match up to the call in "eigenvectors"
  void testEigensystem();

  // the vertex ordering is:
  //
  //   3                  2
  //      o ------------ o 
  //      |              |
  //      |              |
  //      |              |
  //      |              |
  //      |              |
  //      |              |
  //      o ------------ o 
  //   0                   1
  VEC2* _vertices[4];
  VEC2 _restPose[4];

  // 2 x 4 matrix, where each column i is a force on vertex i
  MATRIX _forces;
  VECTOR _forceCache;

  // material vertices
  MATRIX _Dm;

  // material model
  MATERIAL* _material;

  // precache some quantities for force computation
  MATRIX  _H[4];
  MATRIX2 _DmH_Inv[4];
  MATRIX _Bmg[4];
  MATRIX _DFDuA[4];

  VECTOR _flatDFDu[4][8];

  // cache condition number for visualization
  Real _conditionNumber;
  Real _maxEig;
  Real _minEig;
};

#endif
