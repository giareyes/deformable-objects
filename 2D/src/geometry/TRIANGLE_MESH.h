#ifndef TRIANGLE_MESH_H
#define TRIANGLE_MESH_H

#include "SETTINGS.h"
//#include "SQUARE.h"
#include "TRIANGLE.h"
#include "MATERIAL.h"
#include <vector>
#include <map>

// TRIANGLE vertex ordering is CLOCKWISE
class TRIANGLE_MESH
{
public:
  TRIANGLE_MESH(const Real poissonsRatio = 0.3, const Real youngsModulus = 1e6);
  ~TRIANGLE_MESH();

  void buildBeam();
  void buildSingleTest();
  void buildSquashedSquare();
  void buildStretchTest(const int xSquares = 10, const int ySquares = 10);
  void buildSVDTest(const int xSquares = 10, const int ySquares = 10);
  void buildPartialStretchTest(const int xSquares = 10, const int ySquares = 10);
  void buildPullTest(const int xSquares = 10, const int ySquares = 10);
  void buildForceTest(const int xSquares = 10, const int ySquares = 10);
  void buildFloatingTest(const int xSquares = 10, const int ySquares = 10);
  void buildHangTest(const int xSquares = 10, const int ySquares = 10);
  void buildDegenerateHangTest(const int xSquares = 10, const int ySquares = 10, const Real degeneracy = 0.0);
  void buildDegenerateVertexPullTest(const int xSquares = 10, const int ySquares = 10, const Real degeneracy = 0.0);
  void buildDegenerateVertexPullTestBottom(const int xSquares = 10, const int ySquares = 10, const Real degeneracy = 0.0);
  void buildDegenerateAsymmetricVertexPullTest(const int xSquares = 10, const int ySquares = 10, const Real degeneracy = 0.0);
  void buildDegenerateColumnPullTest(const int xSquares = 10, const int ySquares = 10, const Real degeneracy = 0.0);
  void draw() { draw(VEC4(254.0 / 255.0, 240.0 / 255, 217.0 / 255.0, 1.0)); };
  //void draw() { draw(VEC4(0.0, 1.0, 0.0, 1.0)); };
  void drawOutline(const Real lineWidth = 0.0025);
  void draw(const VEC4& color);
  void drawCheckerboard();
  void drawSprings();
  void drawPullSlabs() { drawLeftSlab(); drawRightSlab(); };
  void drawLeftSlab();
  void drawRightSlab();
  void drawForces();
  void drawConstrainedVertices();
  void drawUnconstrainedSurfaceVertices();

  // timestepping routines, will eventually be factored out into
  // their own class
  void stepExplicit();
  void stepImplicit();
  bool stepQuasistatic();
  bool stepQuasistaticWithLineSearch();
  bool stepQuasistaticSandbox();
  bool stepQuasistaticInvertible();
  void addGravity();
  void addBodyForce(const VEC2& bodyForce);
  void reflectMesh();

  // advance the constrained nodes for the stretch test
  void stepStretchTest(const Real stretch);
  void stepSquashTest(const Real squash);
  void stepSpinTest(const Real squash);
  void stepShearTest(const Real stretch);

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
  int pcgIterationsSeen() { return _pcgSeen; };

  // reset material model
  void setMaterial(const std::string& whichMaterial, const Real& poissonsRatio, const Real& youngsModulus);
  void setMaterial(MATERIAL* material);

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
  Real psi(const VECTOR& u);
  Real psiSandbox() const;
  Real psiSandbox(const VECTOR& u);

  // probe the state
  void probeState();

  // normalize condition numbers
  void normalizeConditionNumbers();

  // debugging function -- hand it the solution to 
  // the compression test, see if it is recognized
  void setCompressionSolution(const Real squash);

  Real& springStiffness() { return _springStiffness; };

  // scale the entire mesh
  void scale(const double& scalar);
  
  // scale the current pose in the X direction
  void scaleX(const double& scalar);
  void scaleY(const double& scalar);

  // get a specific triangle
  TRIANGLE& getTriangle(const int index) { return _triangles[index]; };

  const int DOFs() { return _DOFs; };
  const Real conditionNumber() { return _conditionNumber; };
  const Real minEigenvalue() { return _minEig; };
  const Real maxEigenvalue() { return _maxEig; };
  const int pcgIterations() { return _pcgIterations; };
  const Real forceNorm() { return _forceNorm; };

  // get the minimum and maximum interior degrees in the entire mesh
  void minMaxAngles(Real& minDegrees, Real& maxDegrees);
  
  // get the minimum and maximum areas in the entire mesh
  void minMaxAreas(Real& minArea, Real& maxArea);

  // get the worst area quality measure of all the triangles
  Real worstAreaQuality();
  Real worstAngleQuality();
  Real worstKnuppQuality();
  Real worstDmSingularValue();

  // read in a triangle mesh from Shewchuk
  void readShewchukTriangles(const std::string& prefix);

  // build the list of degenerate elements
  void buildDegenerateElementList();

  // return the vertex positions of the unconstrained surface vertices
  VECTOR unconstrainedSurfacePositions();

  // compute which unconstrained vertices are on the surface
  void computeUnconstrainedSurfaceVertices();

  // scramble the vertices to see it untangle
  void scramble();

  // get displacements of all vertices, but ruined vertices stomped
  VECTOR getDisplacementWithoutRuinedVertices();

private:
  // scatter displacement u to the vertices
  void uScatter();
  
  // gather displacement u from the vertices
  void uGather();

  // compute the per-vertex material forces
  void computeMaterialForces();
  void computeMaterialForcesSandbox();
  void computeInvertibleMaterialForces();

  // compute spring constraint forces
  void computeSpringForces();

  // rebuild the vertex-to-index lookup
  void computeVertexToIndexTable();

  // compute the stiffness matrix
  void computeStiffnessMatrix(MATRIX& K);
  void computeStiffnessMatrixSparse(SPARSE_MATRIX& K);
  void computeStiffnessMatrixSparseSandbox(SPARSE_MATRIX& K);
  void computeInvertibleStiffnessMatrixSparse(SPARSE_MATRIX& K);

  // add the spring Jacobians to the stiffness matrix
  void addSpringStiffnesses(SPARSE_MATRIX& K);

  // validate that the stiffness matrix does in fact yield the force
  void validateStiffness(const MATRIX& K, const VECTOR& f);

  // copy the new inversion directions to the old ones
  void updateInversionDirections();

  // return a version with a z coordinate added, just for debugging
  VECTOR padTo3(const VECTOR& v);

  // do a backtracking search in this direction and see if there's a better
  // alpha value
  Real backtrackingLineSearch(const Real original, const VECTOR& direction);
  Real backtrackingLineSearchSandbox(const Real original, const VECTOR& direction);

  // get the minimum/maximum eigenvlue of a sparse matrix
  Real getMaxEig(const SPARSE_MATRIX& K);
  Real getMinEig(const SPARSE_MATRIX& K);
  Real getConditionNumber(const SPARSE_MATRIX& K) { return getMaxEig(K) / getMinEig(K); };

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
  //std::vector<SQUARE> _squares;
  std::vector<TRIANGLE> _triangles;

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
  int _pcgSeen;

  // condition number of first quasistatic iteration
  Real _minEig;
  Real _maxEig;
  Real _conditionNumber;

  // PCG iterations of first quasistatic iteration
  int _pcgIterations;

  // force norm of the first quasistatic iteration
  Real _forceNorm;
  
  // track which elements are degenerate
  std::vector<bool> _isDegenerate;

  // track which unconstrained vertices are on the surface
  std::vector<int> _unconstrainedSurfaceVertices;

  // track wchich vertices were moved to create the degeneracies
  std::vector<int> _ruinedVertices;
};

#endif
