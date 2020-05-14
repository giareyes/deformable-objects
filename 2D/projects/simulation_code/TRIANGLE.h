#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "SETTINGS.h"
#include "MATERIAL.h"
#include "TENSOR4.h"
#include "EXTRAFUNCTIONS.h"
#include <vector>

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

  MATRIX2 computeF() const;

  VEC2 vertexAverage();
  VEC2* vertex(int i) { return _vertices[i]; };
  const VEC2& vertex(int i) const { return *_vertices[i]; };

  MATRIX getConst() { return _constCoef; };
  MATRIX getLinear() { return _linearCoef; };
  TENSOR3 getLinearQuad() { return _quad2; };
  TENSOR3 getQuadCubic() { return _cubic2; };
  TENSOR4 getCubic() { return _cubicCoef; };
  TENSOR4 getQuad() { return _quadraticCoef; };

  vector<MATRIX> computeForceJacobian();

  VECTOR computeForceVector();

  VECTOR precomputedLinearCoef();

  VECTOR precomputedCubicCoef();

  MATRIX precomputedQuadCoef();

  // compute rest area of this triangle
  Real restArea() const;
  Real area() const;

private:
  MATRIX pFpuVectorized();

  std::vector<MATRIX> pFpu();

  VEC2* _vertices[3];
  VEC2 _restPose[3];

  Real _mu;
  Real _lambda;

  MATRIX2 _Dm;
  MATRIX _pfpu;
  MATRIX6 _linearCoef;
  MATRIX6 _constCoef;

  TENSOR4 _quadraticCoef;
  TENSOR3 _cubic2;
  TENSOR3 _quad2;
  TENSOR4 _cubicCoef;

  // material model
  MATERIAL* _material;

  // vector interface for getting and setting positions
  // (only used by unit test)
  VECTOR getDisplacements() const;
  void setDisplacements(const VECTOR& u);
};

#endif
