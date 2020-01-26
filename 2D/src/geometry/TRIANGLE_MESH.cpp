#include "TRIANGLE_MESH.h"
#include "STVK.h"
#include "NEO_HOOKEAN.h"
#include "STABLE_NEO_HOOKEAN.h"
#include "COROTATIONAL.h"
#include "MOONEY_RIVLIN.h"
#include "ARRUDA_BOYCE.h"
#include "VOLUME_NH.h"
#include "VOLUME_MR.h"
#include "DIRICHLET.h"
#include "ARAP.h"
#include "LENGTH_STVK.h"
#include "HYPER_TAN.h"
#include "INVERTIBLE.h"
#include "CARAP.h"
#include "ANISOTROPIC_STVK.h"
#include "ANISOTROPIC_ARAP.h"
#include "COMPOSITE_MATERIAL.h"
#include <iostream>
#include "TIMER.h"
#include "SymEigsSolver.h"
#include "SymEigsShiftSolver.h"
#include "MatOp/SparseSymMatProd.h"
#include "MatOp/SparseSymShiftSolve.h"
#include <float.h>
#include <random>

#define ZEROING_DEGENERATES 0
#define USING_IBEN_ZEROING 0
//#define USING_IBEN_ZEROING 1

#include <GLUT/glut.h>

using namespace std;

TRIANGLE_MESH::TRIANGLE_MESH(const Real poissonsRatio, const Real youngsModulus) : _DOFs(0)
{
  const Real E = youngsModulus;
  const Real nu = poissonsRatio;

  cout << " Young's Modulus: " << youngsModulus << endl;
  cout << " Poisson's Ratio: " << poissonsRatio << endl;

  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const Real mu = E / (2.0 * (1 + nu));
  //_material = new STVK(lambda, mu);
  //_material = new NEO_HOOKEAN(lambda, mu);
  //_material = new CARAP(lambda, mu);
  //_material = new COROTATIONAL(lambda, mu);
  //_material = new ARAP(lambda, mu);
  
  //_material = new ARAP(lambda, mu);
  _material = new STABLE_NEO_HOOKEAN(lambda, mu);
  
  _hasNan = false;
  //_springStiffness = 1000000;
  //_springStiffness = 1000;
  _springStiffness = 10000;
  //_springStiffness = 0;

  _newtonSeen = 0;
  _pcgSeen = 0;
}

TRIANGLE_MESH::~TRIANGLE_MESH() 
{
  delete _material;
}

///////////////////////////////////////////////////////////////////////
// reset material model
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setMaterial(const std::string& whichMaterial, const Real& poissonsRatio, const Real& youngsModulus)
{
  const Real E = youngsModulus;
  const Real nu = poissonsRatio;

  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  //const Real lambda = 0;

  // Some papers use this convention, where incompressibility is then nu = 1.
  // Just to keep things similar-looking between 2D and 3D, we'll use the
  // 3D version.
  //const Real lambda = E * nu / ((1.0 + nu) * (1.0 - nu));
  const Real mu = E / (2.0 * (1 + nu));

  cout << " Using Lame lambda: " << lambda << endl;
  delete _material;
 
  if (strcmp(whichMaterial.c_str(), "STVK") == 0) 
    _material = new STVK(lambda, mu);
  else if (strcmp(whichMaterial.c_str(), "NEO_HOOKEAN") == 0) 
    _material = new NEO_HOOKEAN(lambda, mu);
  else if (strcmp(whichMaterial.c_str(), "COROTATIONAL") == 0) 
    _material = new COROTATIONAL(lambda, mu);
  else if (strcmp(whichMaterial.c_str(), "MOONEY_RIVLIN") == 0) 
    _material = new MOONEY_RIVLIN(lambda, mu);
  else if (strcmp(whichMaterial.c_str(), "ARRUDA_BOYCE") == 0) 
    _material = new ARRUDA_BOYCE(lambda, mu);
  else if (strcmp(whichMaterial.c_str(), "VOLUME_NH") == 0) 
    _material = new VOLUME_NH(lambda, mu);
  else if (strcmp(whichMaterial.c_str(), "ARAP") == 0) 
    _material = new ARAP(lambda, mu);
  else if (strcmp(whichMaterial.c_str(), "INVERTIBLE") == 0)
    _material = new INVERTIBLE(new ARAP(lambda, mu), 0.6);
    //_material = new INVERTIBLE(new STVK(lambda, mu), 0.6);
    //_material = new INVERTIBLE(new NEO_HOOKEAN(lambda, mu), 0.6);
  else if (strcmp(whichMaterial.c_str(), "COMPOSITE") == 0)
  {
    _material = new COMPOSITE_MATERIAL();
    //((COMPOSITE_MATERIAL*)_material)->addMaterial(new VOLUME_NH(lambda, mu));
    //((COMPOSITE_MATERIAL*)_material)->addMaterial(new VOLUME_NH(lambda, 0));
    ((COMPOSITE_MATERIAL*)_material)->addMaterial(new VOLUME_MR(lambda, mu));
    ((COMPOSITE_MATERIAL*)_material)->addMaterial(new DIRICHLET(lambda, mu));
    //((COMPOSITE_MATERIAL*)_material)->addMaterial(new ARAP(lambda, mu));
  }
  else if (strcmp(whichMaterial.c_str(), "HYPER_TAN") == 0)
  {
    _material = new HYPER_TAN(mu, 100 * mu, new ARAP(1,1));
    //_material = new HYPER_TAN(mu, 10 * mu, new ARAP(1,1));
  }
  else
  {
    cout << " MATERIAL TYPE " << whichMaterial.c_str() << " IS UNKNOWN " << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////
// reset material model
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setMaterial(MATERIAL* material)
{
  if (_material)
    delete _material;
  _material = material;

  for (unsigned int x = 0; x < _triangles.size(); x++)
    _triangles[x].setMaterial(_material);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::buildSingleTest()
{
  _vertices.clear();
  _triangles.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  // build the vertices
  VEC2 v0(0,0);
  VEC2 v1(0.5,-sqrt(3.0) / 2.0);
  VEC2 v2(-0.5,-sqrt(3.0) / 2.0);

  _vertices.push_back(v0);
  _vertices.push_back(v1);
  _vertices.push_back(v2);

  _restVertices.push_back(v0);
  _restVertices.push_back(v1);
  _restVertices.push_back(v2);

  _unconstrainedVertices.push_back(0);
  _constrainedVertices.push_back(1);
  _constrainedVertices.push_back(2);

  // build the triangles
  vector<VEC2*> v(4);
  v[0] = &_vertices[0];
  v[1] = &_vertices[1];
  v[2] = &_vertices[2];

  vector<VEC2*> triangle;
  triangle.push_back(v[0]);
  triangle.push_back(v[1]);
  triangle.push_back(v[2]);
  _triangles.push_back(TRIANGLE(_material, triangle));

  // allocate the state vectors
  _DOFs = 2 * (_vertices.size() - _constrainedVertices.size());

  VECTOR zeros(_DOFs);
  zeros.setZero();
  _u          = zeros;
  _uOld       = zeros;
  _f          = zeros;
  _fExternal  = zeros;
  _velocity        = zeros;
  _velocityOld     = zeros;
  _acceleration    = zeros;
  _accelerationOld = zeros;

  // compute the reverse lookup
  computeVertexToIndexTable();
}
/*
{
  _vertices.clear();
  _triangles.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  const int xSquares = 1;
  const int ySquares = 1;

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;

  _springConstraints.clear();

  // build the vertices
  int i = 0;
  for (int y = 0; y < ySquares + 1; y++)
    for (int x = 0; x < xSquares + 1; x++, i++)
    {
      VEC2 vertex(x * xFraction, yStart + y * xFraction);
      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);

      cout << " Vertex: " << vertex.transpose() << endl;

      // vertex 1 is the one not assigned to any mesh
      // vertex 3 is the free vertex on the triangle

      if (i == 0 || i == 2)
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the triangles
  i = 0;
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(4);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];

      vector<VEC2*> upper;
      upper.push_back(v[2]);
      upper.push_back(v[3]);
      upper.push_back(v[0]);
      _triangles.push_back(TRIANGLE(_material, upper));
    }

  // allocate the state vectors
  _DOFs = 2 * (_vertices.size() - _constrainedVertices.size());

  VECTOR zeros(_DOFs);
  zeros.setZero();
  _u          = zeros;
  _uOld       = zeros;
  _f          = zeros;
  _fExternal  = zeros;
  _velocity        = zeros;
  _velocityOld     = zeros;
  _acceleration    = zeros;
  _accelerationOld = zeros;

  // compute the reverse lookup
  computeVertexToIndexTable();
}
*/

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::buildFloatingTest(const int xSquares, const int ySquares)
{
  _vertices.clear();
  _triangles.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;

  _springConstraints.clear();

  // build the vertices
  for (int y = 0; y < ySquares + 1; y++)
    for (int x = 0; x < xSquares + 1; x++)
    {
      VEC2 vertex(x * xFraction, yStart + y * xFraction);
      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);

      //cout << " Vertex: " << vertex.transpose() << endl;

      _unconstrainedVertices.push_back(_vertices.size() - 1);

      /*
      //if (x == xSquares)
      if (x == xSquares)
      {
        pair<VEC2, int> spring;
        spring.first = _vertices.back();

        spring.second = _unconstrainedVertices.size() - 1;
        _springConstraints.push_back(spring);
      }
      */
    }

  // build the triangles
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(4);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];

      vector<VEC2*> upper;
      upper.push_back(v[0]);
      upper.push_back(v[1]);
      upper.push_back(v[2]);
      _triangles.push_back(TRIANGLE(_material, upper));
      
      vector<VEC2*> lower;
      lower.push_back(v[2]);
      lower.push_back(v[3]);
      lower.push_back(v[0]);
      _triangles.push_back(TRIANGLE(_material, lower));

      /*
      if (x == xSquares - 1 && y == 0)
        _lowerRight = _squares.size() - 1;
      if (x == xSquares - 1 && y == ySquares - 1)
        _upperRight = _squares.size() - 1;
      if (x == 0 && y == ySquares - 1)
        _upperLeft = _squares.size() - 1;
        */
    }

  // allocate the state vectors
  _DOFs = 2 * (_vertices.size() - _constrainedVertices.size());

  VECTOR zeros(_DOFs);
  zeros.setZero();
  _u          = zeros;
  _uOld       = zeros;
  _f          = zeros;
  _fExternal  = zeros;
  _velocity        = zeros;
  _velocityOld     = zeros;
  _acceleration    = zeros;
  _accelerationOld = zeros;

  // compute the reverse lookup
  computeVertexToIndexTable();

  this->scale(1.3);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::buildPullTest(const int xSquares, const int ySquares)
{
  _vertices.clear();
  _triangles.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;
  if (yStart < 0) yStart = 0;

  _springConstraints.clear();
  // build the vertices
  for (int y = 0; y < ySquares + 1; y++)
    for (int x = 0; x < xSquares + 1; x++)
    {
      VEC2 vertex(x * xFraction, yStart + y * xFraction);
      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);

      if (x == 0 || x == xSquares)
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the triangles
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(4);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];

      vector<VEC2*> upper;
      upper.push_back(v[0]);
      upper.push_back(v[1]);
      upper.push_back(v[2]);
      _triangles.push_back(TRIANGLE(_material, upper));
      
      vector<VEC2*> lower;
      lower.push_back(v[2]);
      lower.push_back(v[3]);
      lower.push_back(v[0]);
      _triangles.push_back(TRIANGLE(_material, lower));
    }

  // allocate the state vectors
  _DOFs = 2 * (_vertices.size() - _constrainedVertices.size());

  VECTOR zeros(_DOFs);
  zeros.setZero();
  _u          = zeros;
  _uOld       = zeros;
  _f          = zeros;
  _fExternal  = zeros;
  _velocity        = zeros;
  _velocityOld     = zeros;
  _acceleration    = zeros;
  _accelerationOld = zeros;

  // compute the reverse lookup
  computeVertexToIndexTable();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::buildDegenerateColumnPullTest(const int xSquares, const int ySquares, const Real degeneracy)
{
  _vertices.clear();
  _triangles.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;
  if (yStart < 0) yStart = 0;

  std::mt19937 generator(12345);
  std::uniform_real_distribution<double> distribution(-1.0, 1.0);

  _springConstraints.clear();
  // build the square corners 
  for (int y = 0; y < ySquares + 1; y++)
    for (int x = 0; x < xSquares + 1; x++)
    {
      VEC2 vertex(x * xFraction, yStart + y * xFraction);

      bool isDegenerate = false;

      // how much to ruin the conditioning?
      // 0.0 does nothing
      // 1.0 creates all pancakes
      //if (x == (int)(xSquares / 2) && y == 0)
      //if (x == (int)(xSquares / 2) && y < ySquares / 2)
      if (x == (int)(xSquares / 2))
      {
        vertex[0] += xFraction * degeneracy;
        /*
        if (y != 0 && y != ySquares)
          vertex[1] += 0.03 * distribution(generator);
          //vertex[1] += 0.012345 * distribution(generator);
        */
        isDegenerate = true;
      }
      /*
      if (x == (int)(xSquares / 2) && y == 0)
        //vertex[0] += xFraction * 0.5 * degeneracy;
        vertex[0] += xFraction * degeneracy;
      if (x == (int)(xSquares / 2) + 1 && y == 1)
        vertex[0] -= xFraction * 0.5 * degeneracy;
        */

      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);

      if (isDegenerate)
        _ruinedVertices.push_back(_vertices.size() - 1);

      // PULL TEST
      if (x == 0 || x == xSquares)
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the square centers
  const int centerBegin = _vertices.size();
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      VEC2 vertex(x * xFraction + xFraction * 0.5, yStart + y * xFraction + xFraction * 0.5);

      bool isDegenerate = false;

      // how much to ruin the conditioning?
      // 0.0 does nothing
      // 1.0 creates all pancakes
      //if (x == (int)(xSquares / 2) && y == 0)
      if (x == (int)(xSquares / 2))
      //if (x == (int)(xSquares / 2) && y < ySquares / 2)
      {
        vertex[0] += xFraction * 0.5 * degeneracy;
        isDegenerate = true;
        //vertex[0] -= xFraction * 0.5 * degeneracy;
        //vertex[1] += xFraction * 0.5 * degeneracy;
        //vertex[1] -= xFraction * 0.5 * degeneracy;
      }

      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);
      _unconstrainedVertices.push_back(_vertices.size() - 1);
      
      if (isDegenerate)
        _ruinedVertices.push_back(_vertices.size() - 1);
    }

  // build the triangles
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(5);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];
      v[4] = &_vertices[centerBegin + y * xSquares + x];

      vector<VEC2*> pointers;

      // top
      pointers.push_back(v[2]);
      pointers.push_back(v[3]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();

      // left
      pointers.push_back(v[3]);
      pointers.push_back(v[0]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
      
      // bottom
      pointers.push_back(v[0]);
      pointers.push_back(v[1]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
      
      // right
      pointers.push_back(v[1]);
      pointers.push_back(v[2]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
    }

  // allocate the state vectors
  _DOFs = 2 * (_vertices.size() - _constrainedVertices.size());

  VECTOR zeros(_DOFs);
  zeros.setZero();
  _u          = zeros;
  _uOld       = zeros;
  _f          = zeros;
  _fExternal  = zeros;
  _velocity        = zeros;
  _velocityOld     = zeros;
  _acceleration    = zeros;
  _accelerationOld = zeros;

  // compute the reverse lookup
  computeVertexToIndexTable();

  cout << " Total triangles: " << _triangles.size() << endl;
  cout << " Total vertices: " << _vertices.size() << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::buildDegenerateAsymmetricVertexPullTest(const int xSquares, const int ySquares, const Real degeneracy)
{
  _vertices.clear();
  _triangles.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;
  if (yStart < 0) yStart = 0;

  _springConstraints.clear();

  // build the square corners 
  for (int y = 0; y < ySquares + 1; y++)
    for (int x = 0; x < xSquares + 1; x++)
    {
      VEC2 vertex(x * xFraction, yStart + y * xFraction);

      //if (x == (int)(xSquares / 2) - 1 && y == 0)
      //if (x == 2 && y == 0)
      //if (x == (int)(xSquares / 2)  && y == 0)
      if (x == (int)(xSquares / 2))
        vertex[0] += xFraction * degeneracy;

      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);

      // PULL TEST
      if (x == 0 || x == xSquares)
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the triangles
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(4);
      int index = x + y * (xSquares + 1);
      v[0] = &_vertices[index];
      v[1] = &_vertices[index + 1];
      v[2] = &_vertices[index + xSquares + 1];
      v[3] = &_vertices[index + xSquares + 2];

      vector<VEC2*> pointers;

      ///////////////////////////////////////////////////////
      // Upper right to bottom left
      ///////////////////////////////////////////////////////
      // top
      pointers.push_back(v[0]); pointers.push_back(v[1]); pointers.push_back(v[3]);
      // cyclic permutations are valid too:
      //pointers.push_back(v[1]); pointers.push_back(v[3]); pointers.push_back(v[0]);
      //pointers.push_back(v[3]); pointers.push_back(v[0]); pointers.push_back(v[1]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();

      // bottom
      pointers.push_back(v[0]); pointers.push_back(v[3]); pointers.push_back(v[2]);
      // cyclic permutations are valid too:
      //pointers.push_back(v[3]); pointers.push_back(v[2]); pointers.push_back(v[0]);
      //pointers.push_back(v[2]); pointers.push_back(v[0]); pointers.push_back(v[3]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();

      /*
      ///////////////////////////////////////////////////////
      // Upper left to bottom right
      ///////////////////////////////////////////////////////
      // top
      pointers.push_back(v[1]);
      pointers.push_back(v[3]);
      pointers.push_back(v[2]);
      // cyclic permutations are valid too:
      //pointers.push_back(v[3]); pointers.push_back(v[2]); pointers.push_back(v[1]);
      //pointers.push_back(v[2]); pointers.push_back(v[1]); pointers.push_back(v[3]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();

      // bottom
      pointers.push_back(v[0]);
      pointers.push_back(v[1]);
      pointers.push_back(v[2]);
      // cyclic permutations are valid too:
      //pointers.push_back(v[1]); pointers.push_back(v[2]); pointers.push_back(v[0]);
      //pointers.push_back(v[2]); pointers.push_back(v[0]); pointers.push_back(v[1]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
      */
    }

  // allocate the state vectors
  _DOFs = 2 * (_vertices.size() - _constrainedVertices.size());

  VECTOR zeros(_DOFs);
  zeros.setZero();
  _u          = zeros;
  _uOld       = zeros;
  _f          = zeros;
  _fExternal  = zeros;
  _velocity        = zeros;
  _velocityOld     = zeros;
  _acceleration    = zeros;
  _accelerationOld = zeros;

  // compute the reverse lookup
  computeVertexToIndexTable();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::buildDegenerateVertexPullTest(const int xSquares, const int ySquares, const Real degeneracy)
{
  _vertices.clear();
  _triangles.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;
  if (yStart < 0) yStart = 0;

  _springConstraints.clear();
  // build the square corners 
  for (int y = 0; y < ySquares + 1; y++)
    for (int x = 0; x < xSquares + 1; x++)
    {
      VEC2 vertex(x * xFraction, yStart + y * xFraction);
      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);

      // PULL TEST
      if (x == 0 || x == xSquares)
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the square centers
  const int centerBegin = _vertices.size();
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      VEC2 vertex(x * xFraction + xFraction * 0.5, yStart + y * xFraction + xFraction * 0.5);

      bool settingDegenerate = false;
      // how much to ruin the conditioning?
      // 0.0 does nothing
      // 1.0 creates all pancakes
      //if (x == 0 && y == 0)
      if (x == (int)(xSquares / 2) && y == 0)
      //if (x == (int)(xSquares / 2) && y < ySquares / 2 + 1)
      //if (x == (int)(xSquares / 2) && y < ySquares / 2)
      //if (x == (int)(xSquares / 2) && (y == 0 || y == ySquares - 1))
      {
        vertex[0] += xFraction * 0.5 * degeneracy;
        //vertex[0] -= xFraction * 0.5 * degeneracy;
        //vertex[1] += xFraction * 0.5 * degeneracy;
        //vertex[1] -= xFraction * 0.5 * degeneracy;
        //vertex[1] += 0.012345;
        settingDegenerate = true;
      }

      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);
      _unconstrainedVertices.push_back(_vertices.size() - 1);

      if (settingDegenerate)
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " Degenerate vertex index is: " << _vertices.size() - 1 << endl;
        _ruinedVertices.push_back(_vertices.size() - 1);
      }
    }

  // build the triangles
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(5);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];
      v[4] = &_vertices[centerBegin + y * xSquares + x];

      vector<VEC2*> pointers;

      // top
      pointers.push_back(v[2]);
      pointers.push_back(v[3]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();

      // left
      pointers.push_back(v[3]);
      pointers.push_back(v[0]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
      
      // bottom
      pointers.push_back(v[0]);
      pointers.push_back(v[1]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
      
      // right
      pointers.push_back(v[1]);
      pointers.push_back(v[2]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
    }

  // allocate the state vectors
  _DOFs = 2 * (_vertices.size() - _constrainedVertices.size());

  VECTOR zeros(_DOFs);
  zeros.setZero();
  _u          = zeros;
  _uOld       = zeros;
  _f          = zeros;
  _fExternal  = zeros;
  _velocity        = zeros;
  _velocityOld     = zeros;
  _acceleration    = zeros;
  _accelerationOld = zeros;

  // compute the reverse lookup
  computeVertexToIndexTable();
  cout << " Total triangles: " << _triangles.size() << endl;
  cout << " Total vertices: " << _vertices.size() << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::buildDegenerateVertexPullTestBottom(const int xSquares, const int ySquares, const Real degeneracy)
{
  _vertices.clear();
  _triangles.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;
  if (yStart < 0) yStart = 0;

  _springConstraints.clear();
  // build the square corners 
  for (int y = 0; y < ySquares + 1; y++)
    for (int x = 0; x < xSquares + 1; x++)
    {
      VEC2 vertex(x * xFraction, yStart + y * xFraction);

      if (x == (int)(xSquares / 2) && y == 0) 
        vertex[0] += 0.5 * xFraction;
      if (x == (int)(xSquares / 2) + 1 && y == 0)
        vertex[0] -= 0.5 * xFraction;

      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);

      // PULL TEST
      if (x == 0 || x == xSquares)
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the square centers
  const int centerBegin = _vertices.size();
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      VEC2 vertex(x * xFraction + xFraction * 0.5, yStart + y * xFraction + xFraction * 0.5);

      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);
      _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the triangles
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(5);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];
      v[4] = &_vertices[centerBegin + y * xSquares + x];

      vector<VEC2*> pointers;

      // top
      pointers.push_back(v[2]);
      pointers.push_back(v[3]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();

      // left
      pointers.push_back(v[3]);
      pointers.push_back(v[0]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
      
      // bottom
      pointers.push_back(v[0]);
      pointers.push_back(v[1]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
      
      // right
      pointers.push_back(v[1]);
      pointers.push_back(v[2]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
    }

  // allocate the state vectors
  _DOFs = 2 * (_vertices.size() - _constrainedVertices.size());

  VECTOR zeros(_DOFs);
  zeros.setZero();
  _u          = zeros;
  _uOld       = zeros;
  _f          = zeros;
  _fExternal  = zeros;
  _velocity        = zeros;
  _velocityOld     = zeros;
  _acceleration    = zeros;
  _accelerationOld = zeros;

  // compute the reverse lookup
  computeVertexToIndexTable();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::buildHangTest(const int xSquares, const int ySquares)
{
  _vertices.clear();
  _triangles.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;
  if (yStart < 0) yStart = 0;

  _springConstraints.clear();
  // build the square corners 
  for (int y = 0; y < ySquares + 1; y++)
    for (int x = 0; x < xSquares + 1; x++)
    {
      VEC2 vertex(x * xFraction, yStart + y * xFraction);
      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);

      // PULL TEST
      //if (x == 0 || x == xSquares)
      // HANG TEST
      if (x == 0)
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the square centers
  const int centerBegin = _vertices.size();
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      VEC2 vertex(x * xFraction + xFraction * 0.5, yStart + y * xFraction + xFraction * 0.5);

      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);
      _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the triangles
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(5);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];
      v[4] = &_vertices[centerBegin + y * xSquares + x];

      vector<VEC2*> pointers;

      // top
      pointers.push_back(v[2]);
      pointers.push_back(v[3]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();

      // left
      pointers.push_back(v[3]);
      pointers.push_back(v[0]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
      
      // bottom
      pointers.push_back(v[0]);
      pointers.push_back(v[1]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
      
      // right
      pointers.push_back(v[1]);
      pointers.push_back(v[2]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
    }

  // allocate the state vectors
  _DOFs = 2 * (_vertices.size() - _constrainedVertices.size());

  VECTOR zeros(_DOFs);
  zeros.setZero();
  _u          = zeros;
  _uOld       = zeros;
  _f          = zeros;
  _fExternal  = zeros;
  _velocity        = zeros;
  _velocityOld     = zeros;
  _acceleration    = zeros;
  _accelerationOld = zeros;

  // compute the reverse lookup
  computeVertexToIndexTable();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::buildDegenerateHangTest(const int xSquares, const int ySquares, const Real degeneracy)
{
  _vertices.clear();
  _triangles.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;
  if (yStart < 0) yStart = 0;

  _springConstraints.clear();
  // build the square corners 
  for (int y = 0; y < ySquares + 1; y++)
    for (int x = 0; x < xSquares + 1; x++)
    {
      VEC2 vertex(x * xFraction, yStart + y * xFraction);
      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);

      // PULL TEST
      //if (x == 0 || x == xSquares)
      // HANG TEST
      if (x == 0)
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the square centers
  const int centerBegin = _vertices.size();
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      VEC2 vertex(x * xFraction + xFraction * 0.5, yStart + y * xFraction + xFraction * 0.5);

      // how much to ruin the conditioning?
      // 0.0 does nothing
      // 1.0 creates all pancakes
      //if (x == 0 && y == 0)
      if (x == (int)(xSquares / 2) && y == 0)
        vertex[0] += xFraction * 0.5 * degeneracy;

      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);
      _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the triangles
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(5);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];
      v[4] = &_vertices[centerBegin + y * xSquares + x];

      vector<VEC2*> pointers;

      // top
      pointers.push_back(v[2]);
      pointers.push_back(v[3]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();

      // left
      pointers.push_back(v[3]);
      pointers.push_back(v[0]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
      
      // bottom
      pointers.push_back(v[0]);
      pointers.push_back(v[1]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
      
      // right
      pointers.push_back(v[1]);
      pointers.push_back(v[2]);
      pointers.push_back(v[4]);
      _triangles.push_back(TRIANGLE(_material, pointers));
      pointers.clear();
    }

  // allocate the state vectors
  _DOFs = 2 * (_vertices.size() - _constrainedVertices.size());

  VECTOR zeros(_DOFs);
  zeros.setZero();
  _u          = zeros;
  _uOld       = zeros;
  _f          = zeros;
  _fExternal  = zeros;
  _velocity        = zeros;
  _velocityOld     = zeros;
  _acceleration    = zeros;
  _accelerationOld = zeros;

  // compute the reverse lookup
  computeVertexToIndexTable();
}

///////////////////////////////////////////////////////////////////////
// rebuild the vertex-to-index lookup
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeVertexToIndexTable()
{
  _vertexToIndex.clear();
  for (unsigned int x = 0; x < _unconstrainedVertices.size(); x++)
  {
    VEC2* vertex = &_vertices[_unconstrainedVertices[x]];
    _vertexToIndex[vertex] = 2 * x;
  }
}

///////////////////////////////////////////////////////////////////////
// scramble the vertices to see it untangle
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::scramble()
{
  std::mt19937 generator(123);
  std::uniform_real_distribution<double> rando(-2.0, 2.0);

  for (unsigned int x = 0; x < _unconstrainedVertices.size(); x++)
  {
    VEC2& vertex = _vertices[_unconstrainedVertices[x]];
    
    vertex[0] += rando(generator);
    vertex[1] += rando(generator);
  }
  uGather();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::drawSprings()
{
  // draw the springs
  for (int x = 0; x < _springConstraints.size(); x++)
  {
    int meshIndex = _unconstrainedVertices[_springConstraints[x].second];
    VEC2 meshVertex = _vertices[meshIndex];
    glLineWidth(2.0);
    glBegin(GL_LINES);
      glColor4f(0,0,0,1);
      glVertex2f(_springConstraints[x].first[0], _springConstraints[x].first[1]);
      glVertex2f(meshVertex[0], meshVertex[1]);
    glEnd();
    glBegin(GL_POINTS);
      glColor4f(0,0,1,1);
      glVertex2f(_springConstraints[x].first[0], _springConstraints[x].first[1]);
      glColor4f(0,1,0,1);
      glVertex2f(meshVertex[0], meshVertex[1]);
    glEnd();
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::drawUnconstrainedSurfaceVertices()
{
  glColor4f(0,1,0,1);
  glBegin(GL_POINTS);
  for (unsigned int x = 0; x < _unconstrainedSurfaceVertices.size(); x++)
  {
    int index = _unconstrainedSurfaceVertices[x];
    glVertex2f(_vertices[index][0], _vertices[index][1]);
  }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::drawForces()
{
  // draw the forces
  glBegin(GL_LINES);
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    for (int y = 0; y < 3; y++)
    {
      // make sure it's not constrained
      VEC2* vPointer = _triangles[x].vertex(y);
      if (_vertexToIndex.find(vPointer) == _vertexToIndex.end())
        continue;

      VEC2 begin = *vPointer;
      glColor4f(1,0,0,1);
      glVertex2f(begin[0], begin[1]);

      // populate the vector
      int index = _vertexToIndex[vPointer];

      VEC2 distance;
      distance[0] = _f[index];
      distance[1] = _f[index + 1];
      distance.normalize();

      VEC2 end = begin + distance;
      glVertex2f(end[0], end[1]);
    }
  }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::drawLeftSlab()
{
  /*
  // get the min and max of the mesh bounds
  VEC2 mins, maxs;
  boundingBox(mins, maxs);
  */
  static VEC2 mins;
  static VEC2 maxs;
  static bool firstCall = true;

  if (firstCall)
  {
    boundingBox(mins, maxs);
    firstCall = false;
  }

  //VEC3 wallColor(0,0.5,0.25);
  //VEC3 wallColor(0.025,0.01,0.4);
  VEC3 wallColor(0.7,0.7,0.7);

  glColor4f(wallColor[0], wallColor[1], wallColor[2], 1.0);
  // left slab
  glBegin(GL_QUADS);
    //glColor4f(1.0, 1.0, 1.0, 0.0);
    //glColor4f(1.0, 1.0, 1.0, 1.0);
    glVertex2f(-0.1, mins[1] - 0.1);

    //glColor4f(wallColor[0], wallColor[1], wallColor[2], 1.0);
    glVertex2f(0.0,  mins[1] - 0.1);

    //glColor4f(wallColor[0], wallColor[1], wallColor[2], 1.0);
    glVertex2f(0.0,  maxs[1] + 0.1);

    //glColor4f(1.0, 1.0, 1.0, 0.0);
    //glColor4f(1.0, 1.0, 1.0, 1.0);
    glVertex2f(-0.1, maxs[1] + 0.1);
  glEnd();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::drawRightSlab()
{
  // get the min and max of the mesh bounds
  VEC2 mins, maxs;
  boundingBox(mins, maxs);

  //VEC3 wallColor(0,0.5,0.25);
  //VEC3 wallColor(0.025,0.01,0.4);
  VEC3 wallColor(0.7,0.7,0.7);
  glColor4f(wallColor[0], wallColor[1], wallColor[2], 1.0);

  // right slab
  glBegin(GL_QUADS);
    //glColor4f(wallColor[0], wallColor[1], wallColor[2], 1.0);
    glVertex2f(maxs[0],  mins[1] - 0.1);

    //glColor4f(1.0, 1.0, 1.0, 0.0);
    //glColor4f(1.0, 1.0, 1.0, 1.0);
    glVertex2f(maxs[0] + 0.1,  mins[1] - 0.1);

    //glColor4f(1.0, 1.0, 1.0, 0.0);
    //glColor4f(1.0, 1.0, 1.0, 1.0);
    glVertex2f(maxs[0] + 0.1,  maxs[1] + 0.1);

    //glColor4f(wallColor[0], wallColor[1], wallColor[2], 1.0);
    glVertex2f(maxs[0], maxs[1] + 0.1);
  glEnd();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::drawCheckerboard()
{
  // draw the squares
  VEC4 color1(0.6,0.6,0.6,1);
  VEC4 color0(0.3,0.3,0.3,1);
  VEC4 red(1,0,0,1);
  VEC4 color = color0;
  bool firstColor = true;

  int xCheckers = 10;
  //int xCheckers = 8;

  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    if (x % 4 == 0)
    {
      color = (firstColor) ? color0 : color1;
      firstColor = !firstColor;
    }
    if (x % (4 * xCheckers) == 0)
    {
      color = (firstColor) ? color0 : color1;
      firstColor = !firstColor;
    }

    if (_isDegenerate[x] == true)
      _triangles[x].draw(red);
    else
      _triangles[x].draw(color);
  }

  /*
  // try drawing a quadric instead
  const Real pointSize = 0.01;
  static GLUquadricObj* quadric = gluNewQuadric();
  glColor4f(1,0,0,1);
  for (unsigned int x = 0; x < _constrainedVertices.size(); x++)
  {
    glPushMatrix();
      glTranslatef(_vertices[_constrainedVertices[x]][0], 
                   _vertices[_constrainedVertices[x]][1],
                   0.2);
      gluDisk(quadric, 0.0, pointSize, 20, 1);
    glPopMatrix();
  }
  */
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::drawOutline(const Real lineWidth)
{
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    //_triangles[_triangles.size() - x - 1].draw();
    _triangles[x].drawOutline(lineWidth);
    //_squares[x].drawForces();
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::draw(const VEC4& color)
{
  // draw the squares
  //for (unsigned int x = 1; x < 2; x++)
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    //_triangles[_triangles.size() - x - 1].draw();
    _triangles[x].draw(color);
    //_triangles[x].drawOutline();
    //_squares[x].drawForces();
  }

  /*
  // draw the constrained vertices
  //glPointSize(15.0);
  glPointSize(5.0);
  glColor4f(1,0,0,1);
  glBegin(GL_POINTS);
  for (unsigned int x = 0; x < _constrainedVertices.size(); x++)
    glVertex2f(_vertices[_constrainedVertices[x]][0], 
               _vertices[_constrainedVertices[x]][1]);
  glEnd();
  */

  /*
  // try drawing a quadric instead
  const Real pointSize = 0.01;
  static GLUquadricObj* quadric = gluNewQuadric();
  glColor4f(1,0,0,1);
  for (unsigned int x = 0; x < _constrainedVertices.size(); x++)
  {
    glPushMatrix();
      glTranslatef(_vertices[_constrainedVertices[x]][0], 
                   _vertices[_constrainedVertices[x]][1],
                   0.2);
      gluDisk(quadric, 0.0, pointSize, 20, 1);
    glPopMatrix();
  }
  */

 
 /* 
  // draw the unconstrained vertices
  glColor4f(0,0,1,1);
  glBegin(GL_POINTS);
  for (unsigned int x = 0; x < _unconstrainedVertices.size(); x++)
    glVertex2f(_vertices[_unconstrainedVertices[x]][0], 
               _vertices[_unconstrainedVertices[x]][1]);
  glEnd();
  */
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::drawConstrainedVertices()
{
  // try drawing a quadric instead
  const Real pointSize = 0.01;
  static GLUquadricObj* quadric = gluNewQuadric();
  glColor4f(1,0,0,1);
  for (unsigned int x = 0; x < _constrainedVertices.size(); x++)
  {
    glPushMatrix();
      glTranslatef(_vertices[_constrainedVertices[x]][0], 
                   _vertices[_constrainedVertices[x]][1],
                   0.2);
      gluDisk(quadric, 0.0, pointSize, 20, 1);
    glPopMatrix();
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::uScatter()
{
  for (unsigned int x = 0; x < _unconstrainedVertices.size(); x++)
  {
    int index = _unconstrainedVertices[x];

    VEC2 displace(_u[2 * x], _u[2 * x + 1]);

    if (fabs(displace[0]) > 10.0) displace[0] = 0;
    if (fabs(displace[1]) > 10.0) displace[1] = 0;

    _vertices[index][0] = _restVertices[index][0] + displace[0];
    _vertices[index][1] = _restVertices[index][1] + displace[1];
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::uGather()
{
  for (unsigned int x = 0; x < _unconstrainedVertices.size(); x++)
  {
    int index = _unconstrainedVertices[x];
    //_vertices[index][0] = _restVertices[index][0] + _u[2 * x];
    //_vertices[index][1] = _restVertices[index][1] + _u[2 * x + 1];
    _u[2 * x]     = _vertices[index][0] - _restVertices[index][0];
    _u[2 * x + 1] = _vertices[index][1] - _restVertices[index][1];
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeMaterialForcesSandbox()
{
  TIMER functionTimer(__FUNCTION__);
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    // is this one degenerate?
    VECTOR forces(6);
    if (!_isDegenerate[x])
      forces = _triangles[x].computeForceVectorFast();
    else
#if ZEROING_DEGENERATES
      continue;
#elif USING_IBEN_ZEROING
    forces = _triangles[x].computeForceVectorIben();
#else
    forces = _triangles[x].computeForceVectorForDegenerateElement();
#endif
    //VECTOR forces = _triangles[x].computeForceVectorFVM();
    //VECTOR forces = _triangles[x].computeForceVectorCauchy();
    //VECTOR forces = _triangles[x].computeForceVector();
    //VECTOR forces = _triangles[x].computeForceVectorPK2();

    // scatter the forces
    for (int y = 0; y < 3; y++)
    {
      // make sure it's not constrained
      VEC2* vPointer = _triangles[x].vertex(y);
      if (_vertexToIndex.find(vPointer) == _vertexToIndex.end())
        continue;

      // populate the vector
      int index = _vertexToIndex[vPointer];
      _f[index]     += forces[2 * y];
      _f[index + 1] += forces[2 * y + 1];
    }
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeMaterialForces()
{
  TIMER functionTimer(__FUNCTION__);
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
#if 1
    VECTOR forces = _triangles[x].computeForceVectorFast();
    //VECTOR forces = _triangles[x].computeForceVectorFVM();
    //VECTOR forces = _triangles[x].computeForceVectorCauchy();
    //VECTOR forces = _triangles[x].computeForceVector();
    //VECTOR forces = _triangles[x].computeForceVectorPK2();
#else
    VECTOR forces = _triangles[x].computeForceVectorFast();
    VECTOR dec = _triangles[x].computeForceVectorDEC();
    //VECTOR fvm = _triangles[x].computeForceVectorFVM();
    //VECTOR cauchy = _triangles[x].computeForceVectorCauchy();
    //VECTOR PK2 = _triangles[x].computeForceVectorPK2();
    //VECTOR slow = _triangles[x].computeForceVector();

    //VECTOR test = PK2;
    VECTOR test = dec;
    Real diffNorm = (forces - test).norm();

    if (diffNorm > 1e-4)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " DIFF: " << diffNorm << endl;
      cout << " ground: " << forces.transpose() << endl;
      cout << " test: " << test.transpose() << endl;
      cout << " quotient: " << forces.cwiseQuotient(test).transpose() << endl;
      //exit(0);
    }
    else
    {
      cout << " MATCH ";
    }
#endif

    // scatter the forces
    for (int y = 0; y < 3; y++)
    {
      // make sure it's not constrained
      VEC2* vPointer = _triangles[x].vertex(y);
      if (_vertexToIndex.find(vPointer) == _vertexToIndex.end())
        continue;

      // populate the vector
      int index = _vertexToIndex[vPointer];
      _f[index]     += forces[2 * y];
      _f[index + 1] += forces[2 * y + 1];
    }
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::stepExplicit()
{
  MATRIX M(_DOFs, _DOFs);
  M.setIdentity();

  bool verbose = false;

  const Real dampingAlpha = 0.01;
  MATRIX C = dampingAlpha * M;

  const Real dt = 0.001;
  //const Real dt = 0.0001;
  MATRIX lhs = M + dt * 0.5 * C;

  // set the vertex positions according to u
  uScatter();

  if (verbose)
  {
    cout << " INITIAL " << endl;
    cout << " u:    " << padTo3(_u).transpose() << endl;
  }

  // get the material forces
  _f.setZero();
  computeMaterialForces();
  
  // dt^2 appears from the demoninator of the acceleration finite difference
  VECTOR rhs = M * (2.0 * _u - _uOld) + 
               dt * 0.5 * C * _uOld + 
               dt * dt * (_f + _fExternal);

  if (verbose)
  {
    cout << " f:   " << padTo3(_f).transpose() << endl;
    cout << " rhs: " << padTo3(rhs).transpose() << endl;
  }

  VECTOR uNew = lhs.colPivHouseholderQr().solve(rhs);

  _uOld = _u;
  _u = uNew;
  if (verbose)
  {
    cout << " uNew: " << padTo3(uNew).transpose() << endl;
  }

  uScatter();

  _fExternal.setZero();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addGravity()
{
  for (unsigned int x = 0; x < _fExternal.size() / 2; x++)
  {
    //_fExternal[2 * x + 1] -= 10;
    //_fExternal[2 * x + 1] -= 1;
    _fExternal[2 * x + 1] -= 0.001;
    //_fExternal[2 * x + 1] -= 0.0005;
    //_fExternal[2 * x + 1] -= 0.0001;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addBodyForce(const VEC2& bodyForce)
{
  for (unsigned int x = 0; x < _fExternal.size() / 2; x++)
  {
    _fExternal[2 * x]     += bodyForce[0];
    _fExternal[2 * x + 1] += bodyForce[1];
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::stepImplicit()
{
  TIMER functionTimer(__FUNCTION__);
  const Real dt = 0.01;
  //const Real dt = 0.001;
  //const Real dt = 0.0001;

  bool verbose = false;

  // implicit Newmark
  const Real beta = 0.25;
  const Real gamma = 0.5;
  
  // implicit Euler
  //const Real beta = 0.5;
  //const Real gamma = 1.0;

  const Real alpha1 = 1.0 / (beta * dt * dt);
  const Real alpha2 = 1.0 / (beta * dt);
  const Real alpha3 = (1.0 - 2.0 * beta) / (2.0 * beta);
  const Real alpha4 = gamma / (beta * dt);
  const Real alpha5 = 1.0 - gamma / beta;
  const Real alpha6 = (1.0 - gamma / (2.0 * beta)) * dt;

  //MATRIX M(_DOFs, _DOFs);
  SPARSE_MATRIX M(_DOFs, _DOFs);
  M.setIdentity();

  uScatter();
  _uOld = _u;

  if (verbose)
  {
    cout << " INITIAL " << endl;
    cout << " u:    " << padTo3(_u).transpose() << endl;
  }

  Real residual = 1;
  const int maxNewtonIterations = 10;
  int newtonIterations = 0;
  //Real eps = 1e-4;
  Real eps = 1e-3;

  //MATRIX K(_DOFs, _DOFs);
  SPARSE_MATRIX K(_DOFs, _DOFs);
  SPARSE_MATRIX C;
  while (residual > eps && newtonIterations < maxNewtonIterations)
  {
    cout << " Computing K ... " << flush;
    //computeStiffnessMatrix(K);
    computeStiffnessMatrixSparse(K);
    addSpringStiffnesses(K);
    cout << "done. ";

    // negate K, since we use -force on the RHS
    K *= -1;

    const Real dampingAlpha = 0.001;
    const Real dampingBeta = 0.001;
    //const Real dampingAlpha = 0.01;
    //const Real dampingBeta = 0.01;
    //MATRIX C = dampingAlpha * M + dampingBeta * K;
    if (newtonIterations == 0)
      C = dampingAlpha * M + dampingBeta * K;
    
    // build the right hand side
    _f.setZero();
    computeMaterialForces();
    computeSpringForces();

    VECTOR rhs = M * (alpha1 * (_u - _uOld) - alpha2 * _velocity - alpha3 * _acceleration);
    rhs += C * (alpha4 * (_u - _uOld) + alpha5 * _velocity + alpha6 * _acceleration);
    rhs -= _f;
    rhs -= _springForces;
    rhs -= _fExternal;
    residual = rhs.dot(rhs);
    if (verbose)
    {
      cout << endl << endl;
      cout << " NEWTON STEP: " << newtonIterations << endl;
      cout << " f:   " << padTo3(_f).transpose() << endl;
      cout << " rhs: " << padTo3(rhs).transpose() << endl;
      
      VECTOR massRHS = M * (alpha1 * (_u - _uOld) - alpha2 * _velocity - alpha3 * _acceleration);
      VECTOR dampingRHS = C * (alpha4 * (_u - _uOld) + alpha5 * _velocity + alpha6 * _acceleration);
      cout << " mass rhs:    " << padTo3(massRHS).transpose() << endl;
      cout << " damping rhs: " << padTo3(dampingRHS).transpose() << endl;
    }

    cout << " residual: " << residual << endl;
    cout << " velocity:     " << _velocity.norm() << endl;
    cout << " acceleration: " << _acceleration.norm() << endl;

    // all the matrices needed to be built to compute C anyway,
    // so this is the earliest possible break to save on the solve
    if (residual < eps) 
    {
      cout << " Residual is good enough (eps = " << eps << ")." << endl;
      break;
    }

#if 0
    MATRIX lhs = alpha1 * M + alpha4 * C + K;
    cout << " Solving ... " << flush;
    VECTOR uDelta = lhs.colPivHouseholderQr().solve(rhs);
    cout << "done. " << endl;
#else
    SPARSE_MATRIX lhs = alpha1 * M + alpha4 * C + K;
    cout << " Solving PCG ... " << flush;
    TIMER pcgTimer("PCG");
    //ConjugateGradient<MATRIX> cg;
    //cg.compute(lhs);
    ConjugateGradient<SPARSE_MATRIX> cg;
    cg.compute(lhs);
    VECTOR uDelta= cg.solve(rhs);
    pcgTimer.stop();
    cout << "done. Iterations: " << cg.iterations() << " \t error: " << cg.error() << endl;
#endif
    _u -= uDelta;
    uScatter();

    if (verbose)
    {
      cout << " u:      " << padTo3(_u).transpose() << endl;
      cout << " uDelta: " << padTo3(uDelta).transpose() << endl;
    }
    newtonIterations++;
  }

  _accelerationOld = _acceleration;
  _velocityOld     = _velocity;

  _acceleration = alpha1 * (_u - _uOld) - alpha2 * _velocityOld - alpha3 * _accelerationOld;
  _velocity     = alpha4 * (_u - _uOld) + alpha5 * _velocityOld + alpha6 * _accelerationOld;
  if (verbose)
  {
    cout << " FINAL " << endl;
    cout << " velocity:     " << padTo3(_velocity).transpose() << endl;
    cout << " acceleration: " << padTo3(_acceleration).transpose() << endl;
  }

  uScatter();
  cout << " Newton iterations: " << newtonIterations << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::stepQuasistaticSandbox()
{
  TIMER functionTimer(__FUNCTION__);
  uScatter();
  _uOld = _u;
  cout << " u magnitude: " << _u.norm() << endl;

  Real residual = 1;
  int newtonIterations = 0;
  Real eps = 1e-4;

  int maxNewtonIterations = 200;
  //int maxNewtonIterations = 2;
  Real psiInitial = psiSandbox();
  cout << "Psi initial: " << psiInitial << endl;

  // reset PCG iterations to zero to start the count
  _pcgIterations = 0;

  VECTOR bestVector;
  Real bestSeen = 0;
  
  SPARSE_MATRIX K(_DOFs, _DOFs);
  while (residual > eps && newtonIterations < maxNewtonIterations)
  {
    cout << " Computing K ... " << flush;
    //computeStiffnessMatrixSparse(K);
    computeStiffnessMatrixSparseSandbox(K);
    cout << "done. " << endl;

    // add in any spring forces we may be using
    addSpringStiffnesses(K);

    // negate K, since we use -force on the RHS
    K *= -1;

    /*
#define USING_CONDITION_NUMBER 1
#if USING_CONDITION_NUMBER
    Real minEig = getMinEig(K);
    Real maxEig = getMaxEig(K);
    Real conditionNumber = maxEig / minEig; 
#endif
*/

    // build the right hand side
    _f.setZero();
    computeMaterialForcesSandbox();

    computeSpringForces();

    VECTOR rhs(_DOFs);
    rhs.setZero();
    rhs += _f;
    rhs += _springForces;
    rhs += _fExternal;
    residual = rhs.dot(rhs);

    psiInitial = psiSandbox();
    cout << " residual: " << residual << "\t psi: " << psiInitial << endl;

    if (newtonIterations == 0 || residual < bestSeen)
    {
      bestSeen = residual;
      bestVector = _u;
    }

    // all the matrices needed to be built to compute C anyway,
    // so this is the earliest possible break to save on the solve
    if (residual < eps) break;

#define USING_PCG 1
#if USING_PCG
    cout << " Solving PCG ... " << flush;
    TIMER pcgTimer("PCG");
    ConjugateGradient<SPARSE_MATRIX> cg;
    cg.setMaxIterations(10000);
    cg.setTolerance(1e-8);
    cg.compute(K);
    VECTOR uDelta= cg.solve(rhs);
    pcgTimer.stop();
    cout << "done. Iterations: " << cg.iterations() << " \t error: " << cg.error() << endl;

    _pcgSeen += cg.iterations();

#if USING_CONDITION_NUMBER
    cout << " Condition number: " << conditionNumber << endl;
#endif
    if (newtonIterations == 0)
    {
#if USING_CONDITION_NUMBER
      _minEig = minEig;
      _maxEig = maxEig;
      _conditionNumber = conditionNumber;
#endif
      _forceNorm = rhs.norm();
    }
    if (uDelta.hasNaN())
    {
      cout << " RHS: " << rhs << endl;
      cout << " uDelta: " << uDelta << endl;
      cout << " K: " << endl;
      cout << K << endl;
      exit(0);
    }
    _pcgIterations += cg.iterations();
#else
    cout << " Solving Cholesky ... " << flush;
    TIMER choleskyTimer("Cholesky");
    MATRIX denseK = MATRIX(K);    
    LLT<MATRIX> llt;
    llt.compute(denseK);
    VECTOR uDelta= llt.solve(rhs);
    choleskyTimer.stop();
    
    if (newtonIterations == 0)
    {
#if USING_CONDITION_NUMBER
      _minEig = minEig;
      _maxEig = maxEig;
      _conditionNumber = conditionNumber;
#endif
      _pcgIterations = -1;
      _forceNorm = rhs.norm();
    }
#endif

#if 1
    Real alpha = backtrackingLineSearchSandbox(psiInitial, uDelta);
    cout << " Using alpha: " << alpha << endl;
#else
    Real alpha = 1.0;
#endif

    //_u += uDelta;
    _u += alpha * uDelta;
    uScatter();

    newtonIterations++;
  }
  _newtonSeen += newtonIterations;
  
  // clear out any gravity forces
  _fExternal.setZero();

  // make sure we didn't produce a nan
  for (int x = 0; x < _u.size(); x++)
    if (isnan(_u[x]))
      _hasNan = true;

  if (_hasNan)
    _u = _uOld;

  uScatter();
  cout << " Newton iterations: " << newtonIterations << endl;
  if (newtonIterations == maxNewtonIterations)
  {
    cout << " NEWTON DID NOT CONVERGE " << endl;
    //exit(0);
    //_u = _uOld;
    cout << " Using substep with residual: " << bestSeen << endl;
    cout << " Using vector with magnitude: " << bestVector.norm() << endl;
    _u = bestVector;
    return false;
    //return true;
  }

  cout << " Change in u: " << (_uOld - _u).norm() << endl;

  return true;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::stepQuasistaticWithLineSearch()
{
  TIMER functionTimer(__FUNCTION__);
  uScatter();
  _uOld = _u;

  Real residual = 1;
  int newtonIterations = 0;
  Real eps = 1e-4;
  //Real eps = 1e-8;

#define CAPTURING_WEIRD_EQUILIBRIUM 0

#if CAPTURING_WEIRD_EQUILIBRIUM
  int maxNewtonIterations = 1;
#else
  int maxNewtonIterations = 200;
  //int maxNewtonIterations = 50;
  //int maxNewtonIterations = 3;
#endif
  Real psiInitial = psi();
  cout << "Psi initial: " << psiInitial << endl;

  // reset PCG iterations to zero to start the count
  _pcgIterations = 0;

  VECTOR bestVector;
  Real bestSeen = 0;
  
  //MATRIX K(_DOFs, _DOFs);
  SPARSE_MATRIX K(_DOFs, _DOFs);
  while (residual > eps && newtonIterations < maxNewtonIterations)
  {
    cout << " Computing K ... " << flush;
    computeStiffnessMatrixSparse(K);
    cout << "done. " << endl;

    // add in any spring forces we may be using
    addSpringStiffnesses(K);

    // negate K, since we use -force on the RHS
    K *= -1;

#define USING_CONDITION_NUMBER 1
#if USING_CONDITION_NUMBER
    Real minEig = getMinEig(K);
    Real maxEig = getMaxEig(K);
    Real conditionNumber = maxEig / minEig; 
#endif

    // build the right hand side
    _f.setZero();
    computeMaterialForces();

    computeSpringForces();

    VECTOR rhs(_DOFs);
    rhs.setZero();
    rhs += _f;
    rhs += _springForces;
    rhs += _fExternal;
    residual = rhs.dot(rhs);

    psiInitial = psi();
    cout << " residual: " << residual << "\t psi: " << psiInitial << endl;

    if (newtonIterations == 0 || residual < bestSeen)
    {
      bestSeen = residual;
      bestVector = _u;
    }

    // all the matrices needed to be built to compute C anyway,
    // so this is the earliest possible break to save on the solve
    if (residual < eps) break;

#define USING_PCG 1
#if USING_PCG
    cout << " Solving PCG ... " << flush;
    TIMER pcgTimer("PCG");
    ConjugateGradient<SPARSE_MATRIX> cg;
    cg.setMaxIterations(10000);
    cg.compute(K);
    VECTOR uDelta= cg.solve(rhs);
    pcgTimer.stop();
    cout << "done. Iterations: " << cg.iterations() << " \t error: " << cg.error() << endl;

#if USING_CONDITION_NUMBER
    cout << " Condition number: " << conditionNumber << endl;
#endif
    if (newtonIterations == 0)
    {
#if USING_CONDITION_NUMBER
      _minEig = minEig;
      _maxEig = maxEig;
      _conditionNumber = conditionNumber;
#endif
      _forceNorm = rhs.norm();
    }
    if (uDelta.hasNaN())
    {
      cout << " RHS: " << rhs << endl;
      cout << " uDelta: " << uDelta << endl;
      cout << " K: " << endl;
      cout << K << endl;
      exit(0);
    }
    _pcgIterations += cg.iterations();
#else
    cout << " Solving Cholesky ... " << flush;
    TIMER choleskyTimer("Cholesky");
    MATRIX denseK = MATRIX(K);    
    LLT<MATRIX> llt;
    llt.compute(denseK);
    VECTOR uDelta= llt.solve(rhs);
    choleskyTimer.stop();
    
    if (newtonIterations == 0)
    {
#if USING_CONDITION_NUMBER
      _minEig = minEig;
      _maxEig = maxEig;
      _conditionNumber = conditionNumber;
#endif
      _pcgIterations = -1;
      _forceNorm = rhs.norm();
    }
#endif

    Real alpha = backtrackingLineSearch(psiInitial, uDelta);
    cout << " Using alpha: " << alpha << endl;
    //Real alpha = 1.0;

    //_u += uDelta;
    _u += alpha * uDelta;
    uScatter();

    newtonIterations++;
  }
  _newtonSeen = newtonIterations;
  
  // clear out any gravity forces
  _fExternal.setZero();

  // make sure we didn't produce a nan
  for (int x = 0; x < _u.size(); x++)
    if (isnan(_u[x]))
      _hasNan = true;

  if (_hasNan)
    _u = _uOld;

  uScatter();
  cout << " Newton iterations: " << newtonIterations << endl;
  cout << " Total PCG iterations: " << _pcgIterations << endl;

#if !CAPTURING_WEIRD_EQUILIBRIUM
  if (newtonIterations == maxNewtonIterations)
  {
    cout << " NEWTON DID NOT CONVERGE " << endl;
    //exit(0);
    //_u = _uOld;
    cout << " Using substep with residual: " << bestSeen << endl;
    _u = bestVector;
    return false;
    //return true;
  }
#endif

  cout << " Change in u: " << (_uOld - _u).norm() << endl;

  return true;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::stepQuasistatic()
{
  TIMER functionTimer(__FUNCTION__);
  uScatter();
  _uOld = _u;

  Real residual = 1;
  int newtonIterations = 0;
  Real eps = 1e-4;
  //Real eps = 1e-8;

  //int maxNewtonIterations = 10;
  //int maxNewtonIterations = 50;
  int maxNewtonIterations = 100;
  //int maxNewtonIterations = 3;
  //int maxNewtonIterations = 1000;

  cout << "Psi initial: " << psi() << endl;

  VECTOR bestVector;
  Real bestSeen = 0;
  
  //MATRIX K(_DOFs, _DOFs);
  SPARSE_MATRIX K(_DOFs, _DOFs);
  while (residual > eps && newtonIterations < maxNewtonIterations)
  {
    cout << " Computing K ... " << flush;
    //computeStiffnessMatrix(K);
    computeStiffnessMatrixSparse(K);
    cout << "done. " << endl;

    // add in any spring forces we may be using
    addSpringStiffnesses(K);

    // negate K, since we use -force on the RHS
    K *= -1;

    /*
    // DEBUG: peek at the eigenvalues
    MATRIX fullK = MATRIX(K);
    cout << " K eigenvalues: " << endl;
    cout << fullK.eigenvalues() << endl;
    */

    // build the right hand side
    _f.setZero();
    computeMaterialForces();

    computeSpringForces();

    VECTOR rhs(_DOFs);
    rhs.setZero();
    rhs += _f;
    rhs += _springForces;
    rhs += _fExternal;
    residual = rhs.dot(rhs);

    cout << " residual: " << residual << "\t psi: " << psi() << endl;

    if (newtonIterations == 0 || residual < bestSeen)
    {
      bestSeen = residual;
      bestVector = _u;
    }

    // all the matrices needed to be built to compute C anyway,
    // so this is the earliest possible break to save on the solve
    if (residual < eps) break;

#if 0
    MATRIX lhs = K;

    //cout << " lhs: " << lhs << endl;
    cout << " Solving LU ... " << flush;

    ColPivHouseholderQR<MATRIX> householder(K);
    //VECTOR uDelta= lhs.colPivHouseholderQr().solve(rhs);
    VECTOR uDelta = householder.solve(rhs);
    if (householder.info() != Success)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " LU Solve failed!!! " << endl;
      exit(0);
    }
    cout << "done. " << endl;
#else
    cout << " Solving PCG ... " << flush;
    TIMER pcgTimer("PCG");
    //ConjugateGradient<MATRIX> cg;
    //cg.compute(lhs);
    ConjugateGradient<SPARSE_MATRIX> cg;
    cg.compute(K);
    VECTOR uDelta= cg.solve(rhs);
    pcgTimer.stop();
    cout << "done. Iterations: " << cg.iterations() << " \t error: " << cg.error() << endl;
#endif

    //cout << " delta: " << uDelta << endl;
    if (uDelta.hasNaN())
    {
      cout << " RHS: " << rhs << endl;
      cout << " uDelta: " << uDelta << endl;
      cout << " K: " << endl;
      cout << K << endl;
      exit(0);
    }

    _u += uDelta;
    uScatter();

    newtonIterations++;
  }
  _newtonSeen = newtonIterations;
  
  // clear out any gravity forces
  _fExternal.setZero();

  // make sure we didn't produce a nan
  for (int x = 0; x < _u.size(); x++)
    if (isnan(_u[x]))
      _hasNan = true;

  if (_hasNan)
    _u = _uOld;

  uScatter();
  cout << " Newton iterations: " << newtonIterations << endl;
  if (newtonIterations == maxNewtonIterations)
  {
    cout << " NEWTON DID NOT CONVERGE " << endl;
    //exit(0);
    //_u = _uOld;
    cout << " Using substep with residual: " << bestSeen << endl;
    _u = bestVector;
    return false;
    //return true;
  }

  cout << " Change in u: " << (_uOld - _u).norm() << endl;

  return true;
}

///////////////////////////////////////////////////////////////////////
// compute the stiffness matrix
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeStiffnessMatrix(MATRIX& K)
{
  // allocate and clear
  if (K.rows() != _DOFs || K.cols() != _DOFs)
    K = MATRIX(_DOFs, _DOFs);
  K.setZero();

  for (int x = 0; x < _triangles.size(); x++)
  {
    MATRIX localK = _triangles[x].computeForceJacobian();

    /*
    MATRIX test = _triangles[x].computeForceVectorDEC();
    Real diffNorm = (localK - test).norm();
    if (diffNorm > 1e-4)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " DIFF: " << diffNorm << endl;
      cout << " ground: " << localK  << endl;
      cout << " test: " << test << endl;
      exit(0);
    }
    else
    {
      cout << " MATCH ";
    }
    */

    // get the lookup indices
    int indices[3];
    for (int y = 0; y < 3; y++)
    {
      if (_vertexToIndex.find(_triangles[x].vertex(y)) == _vertexToIndex.end())
        indices[y] = -1;
      else
        indices[y] = _vertexToIndex[_triangles[x].vertex(y)];
    }

    // fill in the diagonal entries
    for (int y = 0; y < 3; y++)
    {
      if (indices[y] == -1) continue;

      int index = indices[y];
      K(index,     index)     += localK(2 * y,     2 * y);
      K(index + 1, index)     += localK(2 * y + 1, 2 * y);
      K(index,     index + 1) += localK(2 * y,     2 * y + 1);
      K(index + 1, index + 1) += localK(2 * y + 1, 2 * y + 1);
    }

    // fill in the off-diagonal entries
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
      {
        if (i == j) continue;
        if (indices[i] == -1) continue;
        if (indices[j] == -1) continue;

        int rowStart = indices[i];
        int colStart = indices[j];

        K(rowStart,     colStart)     += localK(2 * i,     2 * j);
        K(rowStart + 1, colStart)     += localK(2 * i + 1, 2 * j);
        K(rowStart,     colStart + 1) += localK(2 * i,     2 * j + 1);
        K(rowStart + 1, colStart + 1) += localK(2 * i + 1, 2 * j + 1);
      }
  }
}

///////////////////////////////////////////////////////////////////////
// compute the stiffness matrix
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeStiffnessMatrixSparseSandbox(SPARSE_MATRIX& K)
{
  TIMER functionTimer(__FUNCTION__);
  // allocate and clear
  if (K.rows() != _DOFs || K.cols() != _DOFs)
    K = SPARSE_MATRIX(_DOFs, _DOFs);
  K.setZero();

  typedef Eigen::Triplet<Real> TRIPLET;
  vector<TRIPLET> triplets;
  triplets.reserve(_DOFs * 3);

  for (int x = 0; x < _triangles.size(); x++)
  {
    // is this one degenerate?
    MATRIX localK(6,6);
    if (!_isDegenerate[x]) 
      localK = _triangles[x].computeForceJacobianDEC();
    else
#if ZEROING_DEGENERATES
      continue;
#elif USING_IBEN_ZEROING
      localK = _triangles[x].computeForceJacobianIben();
#else
      localK = _triangles[x].computeForceJacobianForDegenerateElement();
#endif

    // get the lookup indices
    int indices[3];
    for (int y = 0; y < 3; y++)
    {
      if (_vertexToIndex.find(_triangles[x].vertex(y)) == _vertexToIndex.end())
        indices[y] = -1;
      else
        indices[y] = _vertexToIndex[_triangles[x].vertex(y)];
    }

    // fill in the diagonal entries
    for (int y = 0; y < 3; y++)
    {
      if (indices[y] == -1) continue;

      int index = indices[y];

      for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
          TRIPLET t(index + i, index + j, localK(2 * y + i, 2 * y + j));
          triplets.push_back(t);
        }
    }

    // fill in the off-diagonal entries
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
      {
        if (i == j) continue;
        if (indices[i] == -1) continue;
        if (indices[j] == -1) continue;

        int rowStart = indices[i];
        int colStart = indices[j];

        for (int m = 0; m < 2; m++)
          for (int n = 0; n < 2; n++)
          {
            TRIPLET t(rowStart + m, colStart + n, localK(2 * i + m, 2 * j + n));
            triplets.push_back(t);
          }
      }
  }
  K.setFromTriplets(triplets.begin(), triplets.end());
}

///////////////////////////////////////////////////////////////////////
// compute the stiffness matrix
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeStiffnessMatrixSparse(SPARSE_MATRIX& K)
{
  TIMER functionTimer(__FUNCTION__);
  // allocate and clear
  if (K.rows() != _DOFs || K.cols() != _DOFs)
    K = SPARSE_MATRIX(_DOFs, _DOFs);
  K.setZero();

  typedef Eigen::Triplet<Real> TRIPLET;
  vector<TRIPLET> triplets;
  triplets.reserve(_DOFs * 3);

  for (int x = 0; x < _triangles.size(); x++)
  {
    //MATRIX localK = _triangles[x].computeForceJacobianFast();
    MATRIX localK = _triangles[x].computeForceJacobianDEC();
    /*
    MATRIX localK = _triangles[x].computeForceJacobianFast();

    MATRIX test = _triangles[x].computeForceJacobianDEC();
    cout << " test: " << endl << test << endl;
    cout << " localK: " << endl << localK<< endl;
    Real diffNorm = (localK - test).norm();
    if (diffNorm > 1e-4)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " DIFF: " << diffNorm << endl;
      cout << " ground: " << localK  << endl;
      cout << " test: " << test << endl;
      exit(0);
    }
    else
    {
      cout << " MATCH ";
    }
    */

    // get the lookup indices
    int indices[3];
    for (int y = 0; y < 3; y++)
    {
      if (_vertexToIndex.find(_triangles[x].vertex(y)) == _vertexToIndex.end())
        indices[y] = -1;
      else
        indices[y] = _vertexToIndex[_triangles[x].vertex(y)];
    }

    // fill in the diagonal entries
    for (int y = 0; y < 3; y++)
    {
      if (indices[y] == -1) continue;

      int index = indices[y];

      for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
          TRIPLET t(index + i, index + j, localK(2 * y + i, 2 * y + j));
          triplets.push_back(t);
        }
    }

    // fill in the off-diagonal entries
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
      {
        if (i == j) continue;
        if (indices[i] == -1) continue;
        if (indices[j] == -1) continue;

        int rowStart = indices[i];
        int colStart = indices[j];

        for (int m = 0; m < 2; m++)
          for (int n = 0; n < 2; n++)
          {
            TRIPLET t(rowStart + m, colStart + n, localK(2 * i + m, 2 * j + n));
            triplets.push_back(t);
          }
      }
  }
  K.setFromTriplets(triplets.begin(), triplets.end());
}

///////////////////////////////////////////////////////////////////////
// get the average vertex position of the mesh
///////////////////////////////////////////////////////////////////////
VEC2 TRIANGLE_MESH::meshMean()
{
  VEC2 mean;
  mean.setZero();
  for (int x = 0; x < _vertices.size(); x++)
    mean += _vertices[x];

  mean *= 1.0 / _vertices.size();
  return mean;
}

///////////////////////////////////////////////////////////////////////
// get the min and max of the mesh bounds
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::boundingBox(VEC2& min, VEC2& max)
{
  min = _vertices[0];
  max = _vertices[0];

  for (int x = 1; x < _vertices.size(); x++)
  {
    min[0] = (_vertices[x][0] < min[0]) ? _vertices[x][0] : min[0];
    min[1] = (_vertices[x][1] < min[1]) ? _vertices[x][1] : min[1];
    max[0] = (_vertices[x][0] > max[0]) ? _vertices[x][0] : max[0];
    max[1] = (_vertices[x][1] > max[1]) ? _vertices[x][1] : max[1];
  }
}

///////////////////////////////////////////////////////////////////////
// set a vertex position
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setVertex(const VEC2& v, int index)
{
  assert(index < _vertices.size());

  _vertices[index] = v;
  uGather();
}

///////////////////////////////////////////////////////////////////////
// get the strain energy of the entire mesh
///////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::psiSandbox(const VECTOR& u)
{
  // backup the current state
  uGather();
  VECTOR uBackup = _u;

  // plug in the candidate state
  Real final = 0;
  _u = u;
  uScatter();
  for (int x = 0; x < _triangles.size(); x++)
  {
    // is this one degenerate?
    if (!_isDegenerate[x])
      final += _triangles[x].psi();
    else
      final += _triangles[x].psiDegenerate();
  }

  // restore the current state
  _u = uBackup;
  uScatter();

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the strain energy of the entire mesh
///////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::psi(const VECTOR& u)
{
  // backup the current state
  uGather();
  VECTOR uBackup = _u;

  // plug in the candidate state
  Real final = 0;
  _u = u;
  uScatter();
  for (int x = 0; x < _triangles.size(); x++)
    final += _triangles[x].psi();

  // restore the current state
  _u = uBackup;
  uScatter();

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the strain energy of the entire mesh
///////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::psi() const
{
  Real final = 0;

  for (int x = 0; x < _triangles.size(); x++)
    final += _triangles[x].psi();

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the strain energy of the entire mesh
///////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::psiSandbox() const
{
  Real final = 0;

  for (int x = 0; x < _triangles.size(); x++)
  {
    // is this one degenerate?
    if (!_isDegenerate[x])
      final += _triangles[x].psi();
    else
      final += _triangles[x].psiDegenerate();
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute spring constraint forces
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeSpringForces()
{
  if (_springForces.size() != _DOFs)
    _springForces.resize(_DOFs);
  _springForces.setZero();

  for (int x = 0; x < _springConstraints.size(); x++)
  {
    int index = _springConstraints[x].second;
    VEC2 diff = _springConstraints[x].first - getUnconstrainedVertex(index);
    diff *= _springStiffness;

    _springForces[2 * index] += diff[0];
    _springForces[2 * index + 1] += diff[1];
  }
}

///////////////////////////////////////////////////////////////////////
// add the spring Jacobians to the stiffness matrix
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addSpringStiffnesses(SPARSE_MATRIX& K)
{
  for (int x = 0; x < _springConstraints.size(); x++)
  {
    int index = 2 * _springConstraints[x].second;
    K.coeffRef(index, index) += -_springStiffness;
    K.coeffRef(index + 1, index + 1) += -_springStiffness;
  }
  K.makeCompressed();
}

///////////////////////////////////////////////////////////////////////
// debugging function -- hand it the solution to 
// the compression test, see if it is recognized
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setCompressionSolution(const Real squash)
{
  for (int x = 0; x < _unconstrainedVertices.size(); x++)
  {
    int index = _unconstrainedVertices[x];
    _vertices[index][0] = _restVertices[index][0] * squash;
  }
  uGather();
}

///////////////////////////////////////////////////////////////////////
// scale the entire mesh
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::scale(const double& scalar)
{
  for (int x = 0; x < _vertices.size(); x++)
    _vertices[x] *= scalar;

  cout << " u before: " << _u.norm() << endl;

  uGather();

  // need to do this, otherwise differencing against old shows up as a velocity
  _uOld = _u;
  cout << " u after: " << _u.norm() << endl;
}

///////////////////////////////////////////////////////////////////////
// scale the entire mesh
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::scaleY(const double& scalar)
{
  for (int x = 0; x < _unconstrainedVertices.size(); x++)
    _vertices[_unconstrainedVertices[x]][1] *= scalar;

  uGather();

  // need to do this, otherwise differencing against old shows up as a velocity
  _uOld = _u;
}

///////////////////////////////////////////////////////////////////////
// scale the entire mesh
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::scaleX(const double& scalar)
{
  for (int x = 0; x < _vertices.size(); x++)
    _vertices[x][0] *= scalar;

  cout << " u before: " << _u.norm() << endl;

  uGather();

  // need to do this, otherwise differencing against old shows up as a velocity
  _uOld = _u;
  cout << " u after: " << _u.norm() << endl;
}

///////////////////////////////////////////////////////////////////////
// return a version with a z coordinate added, just for debugging
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::padTo3(const VECTOR& v)
{
  assert((v.size() % 2) == 0);
  int size3 = v.size() / 2 * 3;

  VECTOR final(size3);
  final.setZero();

  for (int x = 0; x < v.size() / 2; x++)
  {
    final[3 * x] = v[2 * x];
    final[3 * x + 1] = v[2 * x + 1];
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the shear test
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::stepShearTest(const Real shear)
{
  for (int x = 0; x < _constrainedVertices.size(); x++)
  {
    int right = _constrainedVertices[x];
    if (_restVertices[right][0] > 0.5)
      _vertices[right][1] += shear;
  }
}

///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the stretch test
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::stepStretchTest(const Real stretch)
{
  for (int x = 0; x < _constrainedVertices.size(); x++)
  {
    int right = _constrainedVertices[x];
    if (_restVertices[right][0] > 0.5)
      _vertices[right][0] += stretch;
  }
}

///////////////////////////////////////////////////////////////////////
// do a backtracking search in this direction and see if there's a 
// better alpha value
///////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::backtrackingLineSearchSandbox(const Real original, const VECTOR& direction)
{
  Real oldEnergy = original;
  Real newEnergy = oldEnergy;
  int iterations = 0; 
  const int maxIterations = 100;
  cout << " Original energy: " << original << std::endl;

  Real alpha = 1.0;
  Real alphaOld = 1.0;
  //Real alpha = 0.2;
  //Real alphaOld = 0.2;
  while (newEnergy <= oldEnergy && iterations < maxIterations)
  {
    VECTOR uNew = _u + alpha * direction;
    oldEnergy = newEnergy;
    newEnergy = psiSandbox(uNew);
    cout << " Energy: " << newEnergy << "\t Alpha: " << alpha << std::endl;

    if (newEnergy >= oldEnergy)
      return alphaOld;
    else
    {
      alphaOld = alpha;
      alpha *= 0.5;
    }

    iterations++;
  }
  return alpha;
}

///////////////////////////////////////////////////////////////////////
// do a backtracking search in this direction and see if there's a 
// better alpha value
///////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::backtrackingLineSearch(const Real original, const VECTOR& direction)
{
  Real oldEnergy = original;
  Real newEnergy = oldEnergy;
  int iterations = 0; 
  const int maxIterations = 100;
  cout << " Original energy: " << original << std::endl;

  Real alpha = 1.0;
  Real alphaOld = 1.0;
  //Real alpha = 0.2;
  //Real alphaOld = 0.2;
  while (newEnergy <= oldEnergy && iterations < maxIterations)
  {
    VECTOR uNew = _u + alpha * direction;
    oldEnergy = newEnergy;
    newEnergy = psi(uNew);
    cout << " Energy: " << newEnergy << "\t Alpha: " << alpha << std::endl;

    if (newEnergy >= oldEnergy)
      return alphaOld;
    else
    {
      alphaOld = alpha;
      alpha *= 0.5;
    }

    iterations++;
  }
  return alpha;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::getMaxEig(const SPARSE_MATRIX& A)
{
  using namespace Spectra;

  SparseSymMatProd<Real> op(A);

  // Construct eigen solver object, requesting the largest three eigenvalues
  const int howMany = 3;
  //const int sandboxSize = 2 * howMany;
  const int sandboxSize = 10 * howMany;
  SymEigsSolver<Real, LARGEST_MAGN, SparseSymMatProd<Real> > eigs(&op, howMany, sandboxSize);
  //SymEigsSolver<Real, LARGEST_ALGE, SparseSymMatProd<Real> > eigs(&op, howMany, sandboxSize);

  // Initialize and compute
  eigs.init();
  int numberConverged = eigs.compute();

  // Retrieve results
  VECTOR evalues;
  if(eigs.info() == SUCCESSFUL)
    evalues = eigs.eigenvalues();
  else
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " TRIANGLE_MESH::getMaxEig FAILED!!!!! " << endl;
    cout << " Converged values: " << numberConverged << endl;
    cout << " Reason: " << eigs.info() << endl;
    return FLT_MAX;
  }

  return evalues[0];
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::getMinEig(const SPARSE_MATRIX& A)
{
  using namespace Spectra;

  // Construct matrix operation object using the wrapper class DenseGenMatProd
  SparseSymShiftSolve<Real> op(A);

  // Construct eigen solver object, requesting the largest three eigenvalues
  //const int howMany = 3;
  //const int sandboxSize = 10 * howMany;
  const int howMany = 1;
  const int sandboxSize = 100;
  SymEigsShiftSolver<Real, LARGEST_MAGN, SparseSymShiftSolve<Real> > eigs(&op, howMany, sandboxSize, 0.0);

  // Initialize and compute
  eigs.init();
  int numberConverged = eigs.compute();

  // Retrieve results
  VECTOR evalues;
  if(eigs.info() == SUCCESSFUL)
    evalues = eigs.eigenvalues();
  else
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " TRIANGLE_MESH::getMinEig FAILED!!!!! " << endl;
    cout << " Converged values: " << numberConverged << endl;
    cout << " Reason: " << eigs.info() << endl;
    return FLT_MIN;
  }

  return evalues[0];
}

///////////////////////////////////////////////////////////////////////
// get the minimum and maximum interior degress in the entire mesh
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::minMaxAngles(Real& minDegrees, Real& maxDegrees)
{
  for (int x = 0; x < _triangles.size(); x++)
  {
    Real currentMin, currentMax;
    _triangles[x].minMaxAngles(currentMin, currentMax);

    if (x == 0)
    {
      minDegrees = currentMin;
      maxDegrees = currentMax;
    }

    minDegrees = (currentMin < minDegrees) ? currentMin : minDegrees;
    maxDegrees = (currentMax > maxDegrees) ? currentMax : maxDegrees;
  }
}

///////////////////////////////////////////////////////////////////////
// get the minimum and maximum areas in the entire mesh
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::minMaxAreas(Real& minArea, Real& maxArea)
{
  for (int x = 0; x < _triangles.size(); x++)
  {
    Real area = _triangles[x].restArea();

    if (x == 0)
    {
      minArea = area;
      maxArea = area;
    }

    minArea = (area < minArea) ? area : minArea;
    maxArea = (area > maxArea) ? area : maxArea;
  }
}

///////////////////////////////////////////////////////////////////////
// get the worst area quality measure of all the triangles
///////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::worstAreaQuality()
{
  Real worstSeen = _triangles[0].areaQuality();
  for (int x = 1; x < _triangles.size(); x++)
  {
    Real current = _triangles[x].areaQuality();
    worstSeen = (current < worstSeen) ? current : worstSeen;
  }

  return worstSeen;
}

///////////////////////////////////////////////////////////////////////
// get the worst angle quality measure of all the triangles
///////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::worstAngleQuality()
{
  Real worstSeen = _triangles[0].angleQuality();
  for (int x = 1; x < _triangles.size(); x++)
  {
    Real current = _triangles[x].angleQuality();
    worstSeen = (current < worstSeen) ? current : worstSeen;
  }

  return worstSeen;
}

///////////////////////////////////////////////////////////////////////
// get the worst composite quality measure of all the triangles,
// according to Knupp 2003, "Algebraic mesh quality metrics for
// unstructured initial meshes"
///////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::worstKnuppQuality()
{
  Real worstSeen = _triangles[0].knuppQuality();
  for (int x = 1; x < _triangles.size(); x++)
  {
    Real current = _triangles[x].knuppQuality();
    worstSeen = (current < worstSeen) ? current : worstSeen;
  }

  return worstSeen;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::worstDmSingularValue()
{
  Real worstSeen = _triangles[0].minDmSingularValue();
  for (int x = 1; x < _triangles.size(); x++)
  {
    Real current = _triangles[x].minDmSingularValue();
    worstSeen = (current < worstSeen) ? current : worstSeen;
  }

  return worstSeen;
}

///////////////////////////////////////////////////////////////////////
// Read in a *.node file
///////////////////////////////////////////////////////////////////////
static void readNodes(const string& filename, vector<VEC2>& nodes, vector<int>& indices)
{
  // read the nodes file
  FILE* file = NULL;
  file = fopen(filename.c_str(), "r");

  if (file == NULL)
  {
    cout << " File " << filename.c_str() << " does not exist! " << endl;
    exit(0);
  }

  // read in the total number of nodes and other attributes
  int totalNodes = -1;
  int dimension = -1;
  int totalAttributes = -1;
  int totalBoundaryMarkers = -1;
  fscanf(file, "%i %i %i %i", &totalNodes, &dimension, &totalAttributes, &totalBoundaryMarkers);

  cout << " Total nodes: " << totalNodes << endl;
  cout << " Dimension: " << dimension << endl;
  cout << " Attributes: " << totalAttributes << endl;
  cout << " Boundary markers: " << totalBoundaryMarkers << endl;

  assert(dimension == 2);

  for (int x = 0; x < totalNodes; x++)
  {
    // get the vertex position
    int index = -1;
    double position[2];
    fscanf(file, "%i %lf %lf", &index, &(position[0]), &(position[1]));
    cout << " index: " << index << "\t node: " << position[0] << " " << position[1] << endl;

    // store it as a node
    VEC2 node(position[0], position[1]);
    nodes.push_back(node);
    indices.push_back(index);

    // strip off the attributes
    double throwAway;
    for (int y = 0; y < totalAttributes; y++)
      fscanf(file, "%lf", &throwAway);

    // strip off the boundary markers
    for (int y = 0; y < totalBoundaryMarkers; y++)
      fscanf(file, "%lf", &throwAway);
  }
  fclose(file);
}

///////////////////////////////////////////////////////////////////////
// Read in a *.ele file
///////////////////////////////////////////////////////////////////////
static void readElements(const string& filename, const int offset, vector<VEC3I>& triangles)
{
  FILE* file = NULL;
  file = fopen(filename.c_str(), "r");

  if (file == NULL)
  {
    cout << " File " << filename.c_str() << " does not exist! " << endl;
    exit(0);
  }

  int totalTriangles = -1;
  int totalNodesPerTriangle = -1;
  int totalAttributes = -1;
  fscanf(file, "%i %i %i", &totalTriangles, &totalNodesPerTriangle, &totalAttributes);

  cout << " Total triangles: " << totalTriangles << endl;
  cout << " Total nodes in each triangle: " << totalNodesPerTriangle << endl;
  cout << " Total attributes per triangle: " << totalAttributes << endl;
  assert(totalNodesPerTriangle == 3);

  for (int x = 0; x < totalTriangles; x++)
  {
    int triangleIndex = -1;
    int nodeIndices[3];

    fscanf(file, "%i %i %i %i", &triangleIndex, &nodeIndices[0], &nodeIndices[1], &nodeIndices[2]);

    VEC3I triangle(nodeIndices[0], nodeIndices[1], nodeIndices[2]);
    triangle -= VEC3I(offset, offset, offset);
    triangles.push_back(triangle);

    cout << " Triangle " << triangleIndex << ": " << triangle[0] << " " << triangle[1] << " " << triangle[2] << endl;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::reflectMesh()
{
  for (int x = 0; x < _vertices.size(); x++)
    _vertices[x][0] = -_vertices[x][0];

  uGather();
}

///////////////////////////////////////////////////////////////////////
// Read in a Delaunay mesh from Shewchuk's "triangle"
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::readShewchukTriangles(const string& prefix)
{
  // read the nodes file
  string nodeFile = prefix + string(".node");
  vector<VEC2> nodes;
  vector<int> indices;
  readNodes(nodeFile, nodes, indices);

  _restVertices = nodes;
  _vertices = nodes;

  // did the file indices start at 1 or 0? If 1, we need to subtract off one
  // from all the triangles' node indices.
  int offset = (indices[0] == 0) ? 0 : 1;

  // read in the elements file
  string elementFile = prefix + string(".ele");
  vector<VEC3I> triangles;
  readElements(elementFile, offset, triangles);

  for (int x = 0; x < triangles.size(); x++)
  {
    vector<VEC2*> vertices;
    for (int y = 0; y < 3; y++)
      vertices.push_back(&_vertices[triangles[x][y]]);
    _triangles.push_back(TRIANGLE(_material, vertices));
  }

  // constrain some verts
  VEC2 mins, maxs;
  boundingBox(mins, maxs);

  for (int x = 0; x < _vertices.size(); x++)
  {
    if (fabs(_vertices[x][0] - mins[0]) < 1e-4 || 
        fabs(_vertices[x][0] - maxs[0]) < 1e-4)
      _constrainedVertices.push_back(x);
    else
      _unconstrainedVertices.push_back(x);
  }
}

///////////////////////////////////////////////////////////////////////
// build the list of degenerate elements
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::buildDegenerateElementList()
{
  _isDegenerate.clear();
  int totalFound = 0;
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    MATRIX Dm = _triangles[x].Dm();
    JacobiSVD<MATRIX2> svd(Dm, ComputeFullU | ComputeFullV);
    VECTOR sigmas = svd.singularValues();
    bool degenerateFound = false;
    for (int y = 0; y < sigmas.size(); y++)
    {
      //if (sigmas[y] < 1e-8)
      //if (sigmas[y] < 1e-3)
      if (sigmas[y] < 1e-5)
      {
        degenerateFound = true;
        totalFound++;
      }
    }
    _isDegenerate.push_back(degenerateFound);
  }

  cout << " Found " << totalFound << " degenerate elements " << endl;
}

///////////////////////////////////////////////////////////////////////
// compute which unconstrained vertices are on the surface
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeUnconstrainedSurfaceVertices()
{
  // see how many radians are associated with each vertex
  int totalUnconstrained = _unconstrainedVertices.size();

  vector<Real> radiansPerVertex(totalUnconstrained);

  // go through each triangle, compute its interior angles
  for (int x = 0; x < _triangles.size(); x++)
  {
    for (int y = 0; y < 3; y++)
    {
      VEC2* vertex = _triangles[x].vertex(y);

      Real angle = _triangles[x].interiorAngle(y);

      if (_vertexToIndex.find(vertex) == _vertexToIndex.end())
        continue;

      // we want the vertex index, not the index into the 2*N monolithic vector,
      // so divide by 2
      int index = _vertexToIndex[vertex] / 2;

      assert(index < totalUnconstrained);
      radiansPerVertex[index] += angle;
    }
  }

  // see which vertices are less than 2 * M_PI
  _unconstrainedSurfaceVertices.clear();
  for (int x = 0; x < radiansPerVertex.size(); x++)
  {
    if (fabs(radiansPerVertex[x] - 360.0) > 1e-3)
    {
      int index = _unconstrainedVertices[x];
      _unconstrainedSurfaceVertices.push_back(index);
    }
  }

  cout << " Found " << _unconstrainedSurfaceVertices.size() << " unconstrained surface vertices " << endl;
}

///////////////////////////////////////////////////////////////////////
// get displacements of all vertices, but ruined vertices stomped
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::getDisplacementWithoutRuinedVertices()
{
  VECTOR u(_vertices.size() * 2);
  u.setZero();

  map<int, bool> isRuined;
  for (unsigned int x = 0; x < _vertices.size(); x++)
    isRuined[x] = false;
  for (unsigned int x = 0; x < _ruinedVertices.size(); x++)
    isRuined[_ruinedVertices[x]] = true;

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    if (isRuined[x]) continue;
    VEC2 displacement = _vertices[x] - _restVertices[x];
    int x2 = 2 * x;

    u[x2] = displacement[0];
    u[x2 + 1] = displacement[1];
  }

  return u;
}

///////////////////////////////////////////////////////////////////////
// return the vertex positions of the unconstrained surface vertices
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::unconstrainedSurfacePositions()
{
  VECTOR positions(_unconstrainedSurfaceVertices.size() * 2);
  for (unsigned int x = 0; x < _unconstrainedSurfaceVertices.size(); x++)
  {
    int index = _unconstrainedSurfaceVertices[x];
    positions[2 * x]     = _vertices[index][0];
    positions[2 * x + 1] = _vertices[index][1];
  }
  return positions;
}
