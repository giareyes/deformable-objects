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

  // build the deformable model
  void buildBlob(const char* filename, bool reduced, int cols);

  // take a step. The function takes in a boolean to state whether or not we are running
  // a reduced simulation
  bool stepQuasistatic(bool reduced);

//----------------------------------------------------------
  // Euler's equation of motion
  void stepMotion(float dt, const VEC2& outerForce);

  // regular equation of motion
  void setMassMatrix(bool reduction);

  // set U without translation by reading in a file and num of cols
  void basisNoTranslation(const char* filename, int basis_cols);

  // reduced coordinates to displacement
  void qTou();

  // add a force to a single vertex
  void addSingleForce(const VEC2& bodyForce, int vertex);

  // precompute coefficients
  void createCoefs();

//----------------------------------------------------------

  // a few functions for deforming the model slightly
  void stretch2(const Real stretch);
  void stepSquashTest(const Real squash);
  void stepShearTest(const Real stretch);
  void addBodyForce(const VEC2& bodyForce);

  MATERIAL* material() { return _material; };
  VECTOR& fExternal() { return _fExternal; };

  // get a specific triangle
  TRIANGLE& getTriangle(const int index) { return _triangles[index]; };

  // get specific details about the triangle mesh
  const int DOFs() { return _DOFs; };
  const std::vector<TRIANGLE>& triangles() { return _triangles; };
  std::vector<VEC2>& vertices() { return _vertices; };
  std::vector<int>& unconstrainedVertices() { return _unconstrainedVertices; };
  void setDisplacement(int index, float d) { _u[index] = d; };
  VECTOR getDisplacement() {return _u; };
  const std::vector<int>& constrainedVertices() { return _constrainedVertices; };
  vector<WALL>& walls() { return _walls; };

private:
  // scatter displacement u to the vertices
  void uScatter();

  // gather displacement u from the vertices
  void uGather();

  // compute internal force using precomputed coefficients
  void computeMaterialForces();

  // compute internal force manually
  void computeUnprecomputedMaterialForces();

  // rebuild the vertex-to-index lookup
  void computeVertexToIndexTable();

  // compute stiffness matrix using precomputed coefficients
  void computeStiffnessMatrix(MATRIX& K);

  // compute stiffness matrix manually
  void computeUnprecomputedStiffnessMatrix(MATRIX& K);

  // how many degrees of freedom are there?
  int _DOFs;

  VECTOR _u; // the displacement vector
  VECTOR _q; // the reduced displacement vector

  MATRIX _U; // change of basis matrix

  VECTOR _f; // the internal force vector
  VECTOR _fExternal; // the external force vector

  MATRIX _mass; // mass matrix

  VECTOR _velocity; // unreduced velocity
  VECTOR _acceleration; // unreduced acceleration

  VECTOR _rv; // reduced velocity
  VECTOR _ra; // reduced acceleration

  MATRIX _constCoef; // constant matrix for stiffness matrix calculation
  TENSOR3 _quadlinear; // linear coefficient for stiffness matrix calculation
  TENSOR4 _quadraticCoef; // quadratic coefficient for stiffness matrix calculation

  MATRIX _linearCoef; // linear coefficient for internal force calculation
  TENSOR3 _cubicquad; // quadratic coefficient for internal force calculation
  TENSOR4 _cubicCoef; // cubic coefficient for internal force calculation

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
