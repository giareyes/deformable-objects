#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "SETTINGS.h"
#include "MATERIAL.h"
#include <vector>
#include "ANISOTROPIC_DIRICHLET.h"
#include "ANISOTROPIC_ARAP.h"

// TRIANGLE vertex ordering is COUNTER CLOCKWISE
//
//          v1
//          o
//         / \
//        /   \
//    e1 /     \ e0
//      /       \
//     /         \
//    /           \
//   o-------------o
// v2      e2      v0
//
class TRIANGLE 
{
public:
  TRIANGLE(MATERIAL* material, const std::vector<VEC2*>& vertices);

  void draw() { draw(VEC4(254.0 / 255.0, 240.0 / 255, 217.0 / 255.0,1.0)); };
  void draw(const VEC4& color);
  void drawOutline(const Real lineWidth);

  // get the deformation gradient at Xi, the normalized position
  // inside a [-1,-1], [1,1] square
  MATRIX2 computeF() const;

  void scale(const VEC2& scalar);
  void translate(const VEC2& translation);
  void rotate(const Real angle);

  VEC2 vertexAverage();

  MATRIX Dm() const { return _Dm; };

  void applyF();

  void applyToRest(const MATRIX2& transform);

  VEC2* vertex(int i) { return _vertices[i]; };

  // get the strain energy inside the triangle
  Real psi() const;
  Real psiDegenerate() const;
  
  MATRIX computeForceJacobian() { return computeForceJacobianFast(); };
  MATRIX computeForceJacobianFast();
  MATRIX computeForceJacobianDEC();
  MATRIX computeForceJacobianForDegenerateElement();
  MATRIX computeForceJacobianIben();
  VECTOR computeForceVector();
  VECTOR computeForceVectorFast();
  VECTOR computeForceVectorForDegenerateElement();
  VECTOR computeForceVectorIben();
 
  // compute forces the DEC way
  VECTOR computeForceVectorDEC();

  // compute forces the FVM way
  VECTOR computeForceVectorFVM();

  // compute forces the Cauchy strain tensor way
  VECTOR computeForceVectorCauchy();

  // compute forces the PK2 way
  VECTOR computeForceVectorPK2();

  // compute rest area of this triangle
  Real restArea() const;
  Real area() const;

  // what is the interior angle at this vertex?
  Real interiorAngle(const int index) const;
  
  // what are the min and max angles on the interior of this triangle?
  void minMaxAngles(Real& minDegrees, Real& maxDegrees);

  // build the area score from [Knupp 2003]
  Real areaQuality();
  
  // build the angle score from [Knupp 2003]
  Real angleQuality();

  // build the composite score from [Knupp 2003]
  Real knuppQuality() { return angleQuality() * areaQuality(); };

  // get the maximum singular value of the material matrix
  Real minDmSingularValue();

  // build vertices of an equilateral triangle with specific edge lengths
  static std::vector<VEC2> equilateral(const Real edgeLength);
  static void runDegeneracyUnitTests();

  // set the material to something else
  void setMaterial(MATERIAL* material) { _material = material; };

  static Real superGlue;

private:
  MATRIX2 pFpu(const int index);
  MATRIX pFpu();
  MATRIX pFpuDegenerate();
  MATRIX pFpuDegenerate(const MATRIX2& DmInverse);
  MATRIX pFpuDegenerateIsotropic(const MATRIX2& DmInverse);
  VECTOR flatten(const MATRIX2& A);

  // pack the vertices into a matrix
  MATRIX vertexMatrix() const;

  VEC2* _vertices[3];
  VEC2 _restPose[3];

  // 2 x 4 matrix, where each column i is a force on vertex i
  MATRIX _forces;

  // material vertices
  MATRIX _Dm;

  // material model
  MATERIAL* _material;

  Real computeIsotropicDegeneratePsi() const;
  Real computeAnisotropicDegeneratePsi() const;

  VECTOR computeIsotropicDegenerateForce();
  VECTOR computeAnisotropicDegenerateForce();
  
  MATRIX computeIsotropicDegenerateForceJacobian();
  MATRIX computeAnisotropicDegenerateForceJacobian();

  void computeDegenerateQuantities(MATRIX2& DmIsotropic,
                                   MATRIX2& DsIsotropic,
                                   MATRIX2& DsAnisotropic,
                                   VEC2& cBar, 
                                   Real& phantomRestArea,
                                   MATERIAL* materialPointer,
                                   MATRIX6& permutation, 
                                   bool debug = false) const;

  // find out which edge is the least degenerate edge
  int findLeastDegenerateEdge() const;

  // analyze a degenerate Dm matrix
  void analyzeDm(const MATRIX2& Dm) const;
  
  // analyze a degenerate Ds matrix
  void analyzeDs(const MATRIX2& Ds) const;

  // terms for taking the full derivative of ansiotropic Ds
  Real cTcPartial(const int index, const int goodEdge);
  MATRIX2 CPartial(const int index, const int goodEdge);
  MATRIX2 DsPartial(const int index);
  MATRIX pFpuDegenerateComplete(const MATRIX2& DmInverse);

  // unit test for a degenerate triangle
  bool testDegenerateDFDu();

  // vector interface for getting and setting positions
  // (only used by unit test)
  VECTOR getDisplacements() const;
  void setDisplacements(const VECTOR& u);

  //
  /*
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
  */
};

#endif
