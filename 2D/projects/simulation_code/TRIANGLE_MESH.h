#ifndef TRIANGLE_MESH_H
#define TRIANGLE_MESH_H

#include "SETTINGS.h"
#include "TRIANGLE.h"
#include "MATERIAL.h"
#include "WALL.h"
#include <vector>
#include <map>

// TRIANGLE vertex ordering is CLOCKWISE
class TRIANGLE_MESH
{
public:
  TRIANGLE_MESH(const Real poissonsRatio = 0.3, const Real youngsModulus = 1e6);
  ~TRIANGLE_MESH();

  // add walls
  void addWall(const WALL& wall)             { _walls.push_back(wall); };
  // build the different kinds of tests
  void buildBlob(const Real xPos, int sceneNum);

  bool stepQuasistatic();

//----------------------------------------------------------
  // Euler's equation of motion
  void stepMotion(float dt, const VEC2& outerForce);

  // regular equation of motion
  void setMassMatrix();

  // set U
  void setBasisReduction();

  void basisNoTranslation();

  void uToq();

  void qTou();

  void addSingleForce(const VEC2& bodyForce, int vertex);
  // your average joe collision detection
  void checkCollision();

  // a weird but fun collision detector
  void wackyCollision();

  // precompute coefficients
  void createCoefs();

//----------------------------------------------------------

  // advance the constrained nodes for the stretch test
  void stepStretchTest(const Real stretch);
  void stretch2(const Real stretch);
  void stepSquashTest(const Real squash);
  void stepShearTest(const Real stretch);
  void addBodyForce(const VEC2& bodyForce);

  MATERIAL* material() { return _material; };
  VECTOR& fExternal() { return _fExternal; };

  // get a specific triangle
  TRIANGLE& getTriangle(const int index) { return _triangles[index]; };

  const int DOFs() { return _DOFs; };
  const std::vector<TRIANGLE>& triangles() { return _triangles; };
  std::vector<VEC2>& vertices() { return _vertices; };
  std::vector<int>& unconstrainedVertices() { return _unconstrainedVertices; };
  void setDisplacement(int index, float d) { _u[index] = d; };
  const std::vector<int>& constrainedVertices() { return _constrainedVertices; };
  vector<WALL>& walls() { return _walls; };

private:
  // scatter displacement u to the vertices
  void uScatter();

  // gather displacement u from the vertices
  void uGather();

  void computeMaterialForces();

  // rebuild the vertex-to-index lookup
  void computeVertexToIndexTable();

  void computeStiffnessMatrix(MATRIX& K);

  // how many degrees of freedom are there?
  int _DOFs;

  // the displacement vector
  VECTOR _u;
  VECTOR _q;

  // change of basis matrix
  MATRIX _U;

  // the force vector
  VECTOR _f;
  VECTOR _fExternal;

  // mass matrix
  MATRIX _mass;

  // velocity and acceleration
  VECTOR _velocity;
  VECTOR _acceleration;

  // reduced velocity and accleration
  VECTOR _rv;
  VECTOR _ra;

  // coefficients for precomputation
  MATRIX _linearCoef;
  MATRIX _constCoef;

  TENSOR4 _quadraticCoef;
  TENSOR3 _quadlinear;
  TENSOR3 _cubicquad;
  TENSOR4 _cubicCoef;

  // the geometry
  std::vector<VEC2>   _vertices;
  std::vector<VEC2>   _restVertices;
  std::vector<TRIANGLE> _triangles;

  // these vertices are constrained
  std::vector<int> _constrainedVertices;

  // these vertices are unconstrained
  std::vector<int> _unconstrainedVertices;

  // back-index from a vertex to its position in _u
  std::map<VEC2*, int> _vertexToIndex;
  std::map<VEC2*, int> _allVertsToIndex;

  // the material
  MATERIAL* _material;

  //walls
  vector<WALL>     _walls;
  int             _flag;
};

#endif
