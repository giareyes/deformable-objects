#ifndef SQUARE_MESH_H
#define SQUARE_MESH_H

#include "SETTINGS.h"
#include "SQUARE.h"
#include "MATERIAL.h"
#include <vector>
#include <map>

class SQUARE_MESH
{
public:
  SQUARE_MESH(const Real poissonsRatio = 0.3, const Real youngsModulus = 1e6);
  ~SQUARE_MESH();

  void buildBeam();
  void buildSquashedSquare();
  void buildStretchTest(const int xSquares = 10, const int ySquares = 10);
  void buildSVDTest(const int xSquares = 10, const int ySquares = 10);
  void buildPartialStretchTest(const int xSquares = 10, const int ySquares = 10);
  void buildGravityTest(const int xSquares = 10, const int ySquares = 10);
  void buildPullTest(const int xSquares = 10, const int ySquares = 10);
  void buildForceTest(const int xSquares = 10, const int ySquares = 10);
  void buildFloatingTest(const int xSquares = 10, const int ySquares = 10);
  void draw();
  void drawSprings();
  void drawInversionDirections();

  // timestepping routines, will eventually be factored out into
  // their own class
  void stepExplicit();
  void stepImplicit();
  bool stepQuasistatic();
  bool stepQuasistaticInvertible();
  void addGravity();

  // advance the constrained nodes for the stretch test
  void stepStretchTest(const Real stretch);
  void stepSquashTest(const Real squash);
  void stepSpinTest(const Real squash);
  void stepShearTest(const Real stretch);
  void stepCycleTest(const Real time);

  // advance the spring anchors for various tests
  void stepStretchForceTest(const Real stretch);
  void stepSquashForceTest(const Real stretch);
  void stepShearForceTest(const Real stretch);

  // get the average vertex position of the mesh
  VEC2 meshMean();

  // get the area of the entire mesh
  Real meshArea() const;

  // get the area of the smallest square
  Real smallestSquareArea() const;
  Real largestSquareArea() const;

  // original area of the mesh
  const Real restArea() const { return _restArea; };

  // ratio of current area to original area
  const Real areaRatio() const { return meshArea() / _restArea; };

  // get the min and max of the mesh bounds
  void boundingBox(VEC2& min, VEC2& max);

  //const MATERIAL* material() { return _material; };
  MATERIAL* material() { return _material; };
  const bool& hasNaN() const { return _hasNan; };
  VECTOR& fExternal() { return _fExternal; };
  int newtonIterationsSeen() { return _newtonSeen; };

  // reset material model
  void setMaterial(const std::string& whichMaterial, const Real& poissonsRatio, const Real& youngsModulus);

  // set a vertex position
  void setVertex(const VEC2& v, int index);
  VEC2 getVertex(int index) { return _vertices[index]; };
  VEC2 getUnconstrainedVertex(int index) { return _vertices[_unconstrainedVertices[index]]; };
  void setUnconstrainedVertex(VEC2& v, int index) { _vertices[_unconstrainedVertices[index]] = v; };
  int totalUnconstrainedVertices() const { return _unconstrainedVertices.size(); };

  VEC2 getConstrainedVertex(int index) { return _vertices[_constrainedVertices[index]]; };
  void setConstrainedVertex(VEC2& v, int index) { _vertices[_constrainedVertices[index]] = v; };
  int totalConstrainedVertices() const { return _constrainedVertices.size(); };

  // get the strain energy of the entire mesh
  Real psi() const;

  // probe the state
  void probeState();

  // normalize condition numbers
  void normalizeConditionNumbers();

  // debugging function -- hand it the solution to 
  // the compression test, see if it is recognized
  void setCompressionSolution(const Real squash);

  // run some SVD-related tests on the first square
  void svdTest() { _squares[0].svdTest(); };

  Real& springStiffness() { return _springStiffness; };

private:
  // scatter displacement u to the vertices
  void uScatter();
  
  // gather displacement u from the vertices
  void uGather();

  // compute the per-vertex material forces
  void computeMaterialForces();
  void computeInvertibleMaterialForces();

  // compute spring constraint forces
  void computeSpringForces();

  // rebuild the vertex-to-index lookup
  void computeVertexToIndexTable();

  // compute the stiffness matrix
  void computeStiffnessMatrix(MATRIX& K);
  void computeStiffnessMatrixSparse(SPARSE_MATRIX& K);
  void computeInvertibleStiffnessMatrixSparse(SPARSE_MATRIX& K);

  // add the spring Jacobians to the stiffness matrix
  void addSpringStiffnesses(SPARSE_MATRIX& K);

  // validate that the stiffness matrix does in fact yield the force
  void validateStiffness(const MATRIX& K, const VECTOR& f);

  // copy the new inversion directions to the old ones
  void updateInversionDirections();

  // how many degrees of freedom are there?
  int _DOFs;

  // the displacement vector
  VECTOR _u;
  VECTOR _uOld;

  // the force vector
  VECTOR _f;
  VECTOR _fExternal;
  VECTOR _velocity;
  VECTOR _velocityOld;
  VECTOR _acceleration;
  VECTOR _accelerationOld;
  VECTOR _springForces;

  // the geometry
  std::vector<VEC2>   _vertices;
  std::vector<VEC2>   _restVertices;
  std::vector<SQUARE> _squares;

  // spring constraints
  // the pair ordering is <Anchor position, unconstrained vertex index>
  std::vector<std::pair<VEC2, int> > _springConstraints;
  Real _springStiffness;

  // sim produced a NaN
  bool _hasNan;

  // areas, so we can check for preservation
  Real _restArea;

  // these vertices are constrained
  std::vector<int> _constrainedVertices;
  
  // these vertices are unconstrained
  std::vector<int> _unconstrainedVertices;

  // back-index from a vertex to its monolithic vector position
  std::map<VEC2*, int> _vertexToIndex;

  // these vertices have forces applied to them
  std::vector<int> _forceVertices;

  // the material
  MATERIAL* _material;

  // specific squares to track
  int _lowerRight;
  int _upperRight;
  int _upperLeft;

  // previous number of newton iterations
  int _newtonSeen;
};

#endif
