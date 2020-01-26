#include "SQUARE_MESH.h"
#include "STVK.h"
#include "NEO_HOOKEAN.h"
#include "STABLE_NEO_HOOKEAN.h"
#include "SCALED_NEO_HOOKEAN.h"
#include "STABLE_ARRUDA_BOYCE.h"
#include "SCALED_ARRUDA_BOYCE.h"
#include "COROTATIONAL.h"
#include "MOONEY_RIVLIN.h"
#include "ARRUDA_BOYCE.h"
#include "VOLUME_NH.h"
#include "VOLUME_MR.h"
#include "ARAP.h"
#include "LENGTH_STVK.h"
#include "HYPER_TAN.h"
#include "INVERTIBLE.h"
#include "COMPOSITE_MATERIAL.h"
#include "CARAP.h"
#include "DISCRETE_CONFORMAL.h"
#include "SYMMETRIC_DIRICHLET.h"
#include "SYMMETRIC_ARAP.h"
#include <iostream>
#include "TIMER.h"

#include <GLUT/glut.h>

using namespace std;

SQUARE_MESH::SQUARE_MESH(const Real poissonsRatio, const Real youngsModulus) : _DOFs(0)
{
  /*
  //const Real nu = 0.3;
  //const Real nu = 0.325;
  const Real nu = 0.3375;
  //const Real nu = 0.35;
  //const Real nu = 0.4;
  //const Real nu = 0.45;
  //const Real nu = 0.49;
  const Real E  = 1e6;
  //const Real E  = 1e5;
  //const Real E  = 1e4;
  //const Real E  = 1e3;
  */
  const Real E = youngsModulus;
  const Real nu = poissonsRatio;

  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const Real mu = E / (2.0 * (1 + nu));
  //_material = new STVK(lambda, mu);
  //_material = new NEO_HOOKEAN(lambda, mu);
  //_material = new SCALED_NEO_HOOKEAN(lambda, mu);
  //_material = new SCALED_ARRUDA_BOYCE(lambda, mu);
  //_material = new CARAP(lambda, mu);
  //_material = new LENGTH_NH(1.0, 1.0);
  //_material = new SYMMETRIC_DIRICHLET(1.0, 1.0);
  //_material = new SYMMETRIC_ARAP(lambda, mu);
  //_material = new STABLE_ARRUDA_BOYCE(lambda, mu);
  //_material = new COROTATIONAL(lambda, mu);
  //_material = new DISCRETE_CONFORMAL(lambda, mu);
  
  //_material = new STABLE_NEO_HOOKEAN(lambda, mu);
  //_material = new LENGTH_NH(lambda, mu);
  //_material = new SYMMETRIC_DIRICHLET(lambda, mu);
  //_material = new ARAP(lambda, mu);
  //_material = new STVK(lambda, mu);
  //_material = new COROTATIONAL(lambda, mu);
  //_material = new DISCRETE_CONFORMAL(lambda, mu);
  
  MATERIAL* material = new NEO_HOOKEAN(lambda, mu);
  _material = new INVERTIBLE(material, 0.9);

  _hasNan = false;
  //_springStiffness = 1000000;
  //_springStiffness = 1000;
  _springStiffness = 10000;
  //_springStiffness = 0;

  _newtonSeen = 0;
}

SQUARE_MESH::~SQUARE_MESH() 
{
  delete _material;
}

///////////////////////////////////////////////////////////////////////
// reset material model
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::setMaterial(const std::string& whichMaterial, const Real& poissonsRatio, const Real& youngsModulus)
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
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::buildBeam()
{
  _vertices.clear();
  _squares.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  //int ySquares = 1;
  //int xSquares = 1;
  
  int ySquares = 2;
  //int xSquares = 5;
  int xSquares = 10;
  //int xSquares = 2;
  //int xSquares = 1;

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;

  // build the vertices
  for (int y = 0; y < ySquares + 1; y++)
    for (int x = 0; x < xSquares + 1; x++)
    {
      VEC2 vertex(x * xFraction, yStart + y * xFraction);
      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);

      if (x == 0)
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the squares
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(4);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];

      _squares.push_back(SQUARE(_material, v));
    }

  // allocate the state vectors
  _DOFs = 2 * (_vertices.size() - _constrainedVertices.size());

  VECTOR zeros(_DOFs);
  zeros.setZero();
  _u          = zeros;
  _uOld       = zeros;
  _f          = zeros;
  _fExternal  = zeros;
  _springForces    = zeros;
  _velocity        = zeros;
  _velocityOld     = zeros;
  _acceleration    = zeros;
  _accelerationOld = zeros;

  // compute the reverse lookup
  computeVertexToIndexTable();

  // cache the rest area
  _restArea = meshArea();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::buildFloatingTest(const int xSquares, const int ySquares)
{
  _vertices.clear();
  _squares.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;

  // build the vertices
  for (int y = 0; y < ySquares + 1; y++)
    for (int x = 0; x < xSquares + 1; x++)
    {
      VEC2 vertex(x * xFraction, yStart + y * xFraction);
      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);

      /*
      if (x == 0)
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);

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

  // build the squares
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(4);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];

      _squares.push_back(SQUARE(_material, v));

      if (x == xSquares - 1 && y == 0)
        _lowerRight = _squares.size() - 1;
      if (x == xSquares - 1 && y == ySquares - 1)
        _upperRight = _squares.size() - 1;
      if (x == 0 && y == ySquares - 1)
        _upperLeft = _squares.size() - 1;
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
  
  // cache the rest area
  _restArea = meshArea();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::buildForceTest(const int xSquares, const int ySquares)
{
  _vertices.clear();
  _squares.clear();
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

      if (x == 0)
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);

      //if (x == xSquares)
      if (x == xSquares)
      {
        pair<VEC2, int> spring;
        spring.first = _vertices.back();

        spring.second = _unconstrainedVertices.size() - 1;
        _springConstraints.push_back(spring);
      }
    }

  // build the squares
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(4);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];

      _squares.push_back(SQUARE(_material, v));

      if (x == xSquares - 1 && y == 0)
        _lowerRight = _squares.size() - 1;
      if (x == xSquares - 1 && y == ySquares - 1)
        _upperRight = _squares.size() - 1;
      if (x == 0 && y == ySquares - 1)
        _upperLeft = _squares.size() - 1;
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
  
  // cache the rest area
  _restArea = meshArea();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::buildStretchTest(const int xSquares, const int ySquares)
{
  _vertices.clear();
  _squares.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  //int ySquares = 1;
  //int xSquares = 1;

/*  
  int ySquares = 10;
  //int ySquares = 5;
  //int xSquares = 5;
  int xSquares = 10;
  //int xSquares = 2;
  //int xSquares = 1;
 
  //int amp = 8; 
  //int amp = 4; 
  int amp = 1; 
  ySquares *= amp;
  xSquares *= amp;
  */

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;

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

  // build the squares
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(4);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];

      _squares.push_back(SQUARE(_material, v));
      if (x == xSquares - 1 && y == 0)
        _lowerRight = _squares.size() - 1;
      if (x == xSquares - 1 && y == ySquares - 1)
        _upperRight = _squares.size() - 1;
      if (x == 0 && y == ySquares - 1)
        _upperLeft = _squares.size() - 1;
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
  
  // cache the rest area
  _restArea = meshArea();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::buildSVDTest(const int xSquares, const int ySquares)
{
  _vertices.clear();
  _squares.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;

  // build the vertices
  for (int y = 0; y < ySquares + 1; y++)
    for (int x = 0; x < xSquares + 1; x++)
    {
      VEC2 vertex(x * xFraction, yStart + y * xFraction);

      Real theta = M_PI / 4;
      //Real theta = 0;
      MATRIX2 rotate;
      rotate << cos(theta), -sin(theta),
                sin(theta), cos(theta); 
      vertex = rotate * vertex;

      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);

      if (x == 0 || x == xSquares)
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the squares
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(4);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];

      _squares.push_back(SQUARE(_material, v));
      if (x == xSquares - 1 && y == 0)
        _lowerRight = _squares.size() - 1;
      if (x == xSquares - 1 && y == ySquares - 1)
        _upperRight = _squares.size() - 1;
      if (x == 0 && y == ySquares - 1)
        _upperLeft = _squares.size() - 1;
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
  
  // cache the rest area
  _restArea = meshArea();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::buildGravityTest(const int xSquares, const int ySquares)
{
  _vertices.clear();
  _squares.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;

  // build the vertices
  for (int y = 0; y < ySquares + 1; y++)
    for (int x = 0; x < xSquares + 1; x++)
    {
      VEC2 vertex(x * xFraction, yStart + y * xFraction);
      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);

      //if (x == 0 || (x == xSquares && vertex[1] > 0.5 && vertex[1] < 0.4))
      if (y == 0)
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the squares
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(4);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];

      _squares.push_back(SQUARE(_material, v));
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
  
  // cache the rest area
  _restArea = meshArea();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::buildPartialStretchTest(const int xSquares, const int ySquares)
{
  _vertices.clear();
  _squares.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  Real xFraction = 1.0 / xSquares;
  Real yStart = 0.5 - xFraction;

  // build the vertices
  for (int y = 0; y < ySquares + 1; y++)
    for (int x = 0; x < xSquares + 1; x++)
    {
      VEC2 vertex(x * xFraction, yStart + y * xFraction);
      _vertices.push_back(vertex);
      _restVertices.push_back(vertex);

      //if (x == 0 || (x == xSquares && vertex[1] > 0.5 && vertex[1] < 0.4))
      if (x == 0 || (x == xSquares && vertex[1] > 0.8 && vertex[1] < 1.1))
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);
    }

  // build the squares
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(4);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];

      _squares.push_back(SQUARE(_material, v));
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
  
  // cache the rest area
  _restArea = meshArea();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::buildPullTest(const int xSquares, const int ySquares)
{
  _vertices.clear();
  _squares.clear();
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

      //if (x == 0 || x == xSquares)
      //if (x == 0)
      if (y == 0)
        _constrainedVertices.push_back(_vertices.size() - 1);
      else
        _unconstrainedVertices.push_back(_vertices.size() - 1);

      if (y == ySquares && x != 0 && x != xSquares)
      {
        pair<VEC2, int> spring;
        spring.first = _vertices.back();

        spring.first[1] = yStart - 0.25;

        /*
        spring.first[0] -= 0.5;
        spring.first[0] *= 2.0;
        spring.first[0] += 0.5;
        */

        spring.second = _unconstrainedVertices.size() - 1;
        _springConstraints.push_back(spring);
      }
    }

  // build the squares
  for (int y = 0; y < ySquares; y++)
    for (int x = 0; x < xSquares; x++)
    {
      vector<VEC2*> v(4);
      v[0] = &_vertices[y * (xSquares + 1) + x];
      v[1] = &_vertices[y * (xSquares + 1) + x + 1];
      v[2] = &_vertices[(y + 1) * (xSquares + 1) + x + 1];
      v[3] = &_vertices[(y + 1) * (xSquares + 1) + x];

      _squares.push_back(SQUARE(_material, v));
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
  
  // cache the rest area
  _restArea = meshArea();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::buildSquashedSquare()
{
  _vertices.clear();
  _squares.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  _vertices = vector<VEC2>(4);
  _vertices[0] = VEC2(-1, -1);
  _vertices[1] = VEC2(1, -1);
  _vertices[2] = VEC2(1, 1);
  _vertices[3] = VEC2(-1, 1);

  for (int x = 0; x < _vertices.size(); x++)
  {
    _restVertices.push_back(_vertices[x]);
    _unconstrainedVertices.push_back(x);
  }

  // build the square
  vector<VEC2*> v(4);
  v[0] = &_vertices[0];
  v[1] = &_vertices[1];
  v[2] = &_vertices[2];
  v[3] = &_vertices[3];

  _squares.push_back(SQUARE(_material, v));

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

  // squash the square
  _squares[0].scale(VEC2(1.3, 1));
  _squares[0].scale(VEC2(0.5, 1.5));
  _squares[0].rotate(M_PI / 8);

  _u = _squares[0].getDisplacement();
  //(*_vertices[0])[0] += 0.1;
}

///////////////////////////////////////////////////////////////////////
// rebuild the vertex-to-index lookup
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::computeVertexToIndexTable()
{
  _vertexToIndex.clear();
  for (unsigned int x = 0; x < _unconstrainedVertices.size(); x++)
  {
    VEC2* vertex = &_vertices[_unconstrainedVertices[x]];
    _vertexToIndex[vertex] = 2 * x;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::drawInversionDirections()
{
  for (int x = 0; x < _squares.size(); x++)
    _squares[x].drawInversionDirection();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::drawSprings()
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
// normalize condition numbers
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::normalizeConditionNumbers()
{
  static int frames = 0;
  static Real minFound = 0;
  static Real normalize = 0;

  // only do this the first time
  //if (frames <= 10)
  {
    Real maxFound = _squares[0].conditionNumber();
    minFound = _squares[0].conditionNumber();
    for (unsigned int x = 1; x < _squares.size(); x++)
    {
      maxFound = (maxFound < _squares[x].conditionNumber()) ? _squares[x].conditionNumber() : maxFound;
      minFound = (minFound > _squares[x].conditionNumber()) ? _squares[x].conditionNumber() : minFound;
    }
    Real diff = (maxFound - minFound);
    normalize = (diff > 0) ? 1.0 / diff : 1;

    cout << " condition number max: " << maxFound << " min: " << minFound << endl;
  }
  frames++;

  // normalize the condition number
  //cout << " Drawing " << _squares.size() << " squares " << endl;
  //cout << " upper left: " << _upperLeft << endl;
  for (unsigned int x = 0; x < _squares.size(); x++)
  {
    Real number = _squares[x].conditionNumber();
    number -= minFound;
    number *= normalize;
    
    /*
    if (x == _upperLeft)
    {
      cout << " upper left: " << endl;
      cout << " Before: " << _squares[x].conditionNumber() << " after: " << number << endl;
    }
    if (x == _upperRight)
    {
      cout << " upper right: " << endl;
      cout << " Before: " << _squares[x].conditionNumber() << " after: " << number << endl;
    }
    if (x == _lowerRight)
    {
      cout << " lower right: " << endl;
      cout << " Before: " << _squares[x].conditionNumber() << " after: " << number << endl;
    }
    */

    _squares[x].conditionNumber() = number;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::draw()
{
  /*
  static int frames = 0;
  static Real minFound;
  static Real normalize;

  // only do this the first time
  //if (frames <= 10)
  {
    Real maxFound = _squares[0].conditionNumber();
    minFound = _squares[0].conditionNumber();
    for (unsigned int x = 0; x < _squares.size(); x++)
    {
      maxFound = (maxFound < _squares[x].conditionNumber()) ? _squares[x].conditionNumber() : maxFound;
      minFound = (minFound > _squares[x].conditionNumber()) ? _squares[x].conditionNumber() : minFound;
    }
    Real diff = (maxFound - minFound);
    normalize = (diff > 0) ? 1.0 / diff : 1;

    cout << " max: " << maxFound << " min: " << minFound << endl;
  }
  frames++;

  // normalize the condition number
  cout << " Drawing " << _squares.size() << " squares " << endl;
  cout << " upper left: " << _upperLeft << endl;
  for (unsigned int x = 0; x < _squares.size(); x++)
  {
    Real number = _squares[x].conditionNumber();
    number -= minFound;
    number *= normalize;
    
    if (x == _upperLeft)
    {
      cout << " upper left: " << endl;
      cout << " Before: " << _squares[x].conditionNumber() << " after: " << number << endl;
    }
    if (x == _upperRight)
    {
      cout << " upper right: " << endl;
      cout << " Before: " << _squares[x].conditionNumber() << " after: " << number << endl;
    }
    if (x == _lowerRight)
    {
      cout << " lower right: " << endl;
      cout << " Before: " << _squares[x].conditionNumber() << " after: " << number << endl;
    }

    _squares[x].conditionNumber() = number;
  }
  */

  // draw the squares
  for (unsigned int x = 0; x < _squares.size(); x++)
  {
    _squares[x].draw();
    //_squares[x].drawForces();
  }

  // draw the constrained vertices
  //glPointSize(15.0);
  glPointSize(5.0);
  glColor4f(1,0,0,1);
  glBegin(GL_POINTS);
  for (unsigned int x = 0; x < _constrainedVertices.size(); x++)
    glVertex2f(_vertices[_constrainedVertices[x]][0], 
               _vertices[_constrainedVertices[x]][1]);
  glEnd();
 
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
void SQUARE_MESH::uScatter()
{
  for (unsigned int x = 0; x < _unconstrainedVertices.size(); x++)
  {
    int index = _unconstrainedVertices[x];
    _vertices[index][0] = _restVertices[index][0] + _u[2 * x];
    _vertices[index][1] = _restVertices[index][1] + _u[2 * x + 1];
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::uGather()
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
void SQUARE_MESH::computeMaterialForces()
{
  TIMER functionTimer(__FUNCTION__);
  for (unsigned int x = 0; x < _squares.size(); x++)
  {
    //VECTOR forces = _squares[x].computeForceVector();
    VECTOR forces = _squares[x].computeForceVectorFast();

    // scatter the forces
    for (int y = 0; y < 4; y++)
    {
      // make sure it's not constrained
      VEC2* vPointer = _squares[x].vertex(y);
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
void SQUARE_MESH::computeInvertibleMaterialForces()
{
  TIMER functionTimer(__FUNCTION__);
  for (unsigned int x = 0; x < _squares.size(); x++)
  {
    //VECTOR forces = _squares[x].computeForceVector();
    VECTOR forces = _squares[x].computeInvertibleForceVector();

    // scatter the forces
    for (int y = 0; y < 4; y++)
    {
      // make sure it's not constrained
      VEC2* vPointer = _squares[x].vertex(y);
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
void SQUARE_MESH::stepExplicit()
{
  MATRIX M(_DOFs, _DOFs);
  M.setIdentity();

  const Real dampingAlpha = 0.5;
  MATRIX C = dampingAlpha * M;

  const Real dt = 0.01;
  MATRIX lhs = M + dt * 0.5 * C;

  // set the vertex positions according to u
  uScatter();

  // get the material forces
  _f.setZero();
  computeMaterialForces();
  
  // dt^2 appears from the demoninator of the acceleration finite difference
  VECTOR rhs = M * (2.0 * _u - _uOld) + 
               dt * 0.5 * C * _uOld + 
               dt * dt * (_f + _fExternal);

  VECTOR uNew = lhs.colPivHouseholderQr().solve(rhs);

  _uOld = _u;
  _u = uNew;

  uScatter();
  //cout << " displacement magnitude: " << _u.dot(_u) << endl;

  //_fExternal.setZero();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::addGravity()
{
  for (unsigned int x = 0; x < _fExternal.size() / 2; x++)
  {
    _fExternal[2 * x + 1] -= 0.001;
    //_fExternal[2 * x + 1] -= 20;
    //_fExternal[2 * x + 1] -= 1;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::stepImplicit()
{
  const Real dt = 0.01;
  //const Real dt = 0.001;
  //const Real dt = 0.0001;

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

  Real residual = 1;
  const int maxNewtonIterations = 10;
  int newtonIterations = 0;
  Real eps = 1e-4;

  //MATRIX K(_DOFs, _DOFs);
  SPARSE_MATRIX K(_DOFs, _DOFs);
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
    //MATRIX C = dampingAlpha * M + dampingBeta * K;
    SPARSE_MATRIX C = dampingAlpha * M + dampingBeta * K;
    
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

    cout << " residual: " << residual << endl;

    // all the matrices needed to be built to compute C anyway,
    // so this is the earliest possible break to save on the solve
    if (residual < eps) break;

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

    newtonIterations++;
  }

  _accelerationOld = _acceleration;
  _velocityOld     = _velocity;

  _acceleration = alpha1 * (_u - _uOld) - alpha2 * _velocityOld - alpha3 * _accelerationOld;
  _velocity     = alpha4 * (_u - _uOld) + alpha5 * _velocityOld + alpha6 * _accelerationOld;

  uScatter();
  cout << " Newton iterations: " << newtonIterations << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
bool SQUARE_MESH::stepQuasistatic()
{
  TIMER functionTimer(__FUNCTION__);
  uScatter();
  _uOld = _u;

  Real residual = 1;
  int newtonIterations = 0;
  Real eps = 1e-4;
  //Real eps = 1e-8;

  //int maxNewtonIterations = 10;
  int maxNewtonIterations = 50;
  //int maxNewtonIterations = 100;
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

  Real area = meshArea();
  Real ratio = area / _restArea;

  cout << " Current area: " << area 
       << "\t Original area: " << _restArea 
       << "\t Ratio: " << ratio << endl;

  cout << " Change in u: " << (_uOld - _u).norm() << endl;

  return true;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
bool SQUARE_MESH::stepQuasistaticInvertible()
{
  TIMER functionTimer(__FUNCTION__);
  uScatter();
  _uOld = _u;

  Real residual = 1;
  int newtonIterations = 0;
  Real eps = 1e-4;
  //Real eps = 1e-8;

  //int maxNewtonIterations = 10;
  int maxNewtonIterations = 50;
  //int maxNewtonIterations = 100;
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
    //computeStiffnessMatrixSparse(K);
    computeInvertibleStiffnessMatrixSparse(K);
    cout << "done. " << endl;

    // add in any spring forces we may be using
    addSpringStiffnesses(K);

    // negate K, since we use -force on the RHS
    K *= -1;

    /*
    // DEBUG: peek at the eigenvalues
    MATRIX fullK = MATRIX(K);
    cout << "K: " << endl << fullK << endl;
    cout << "K eigenvalues: " << endl;
    cout << fullK.eigenvalues() << endl;
    */

    // build the right hand side
    _f.setZero();
    computeInvertibleMaterialForces();
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

#if 1
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
    VECTOR uDelta = cg.solve(rhs);
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

    updateInversionDirections();

    newtonIterations++;
  }
  _newtonSeen = newtonIterations;

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

  Real area = meshArea();
  Real ratio = area / _restArea;

  cout << " Current area: " << area 
       << "\t Original area: " << _restArea 
       << "\t Ratio: " << ratio << endl;

  cout << " Change in u: " << (_uOld - _u).norm() << endl;

  return true;
}

///////////////////////////////////////////////////////////////////////
// compute the stiffness matrix
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::computeStiffnessMatrix(MATRIX& K)
{
  // allocate and clear
  if (K.rows() != _DOFs || K.cols() != _DOFs)
    K = MATRIX(_DOFs, _DOFs);
  K.setZero();

  for (int x = 0; x < _squares.size(); x++)
  {
    MATRIX localK = _squares[x].computeForceJacobian();

    // get the lookup indices
    int indices[4];
    for (int y = 0; y < 4; y++)
    {
      if (_vertexToIndex.find(_squares[x].vertex(y)) == _vertexToIndex.end())
        indices[y] = -1;
      else
        indices[y] = _vertexToIndex[_squares[x].vertex(y)];
    }

    // fill in the diagonal entries
    for (int y = 0; y < 4; y++)
    {
      if (indices[y] == -1) continue;

      int index = indices[y];
      K(index,     index)     += localK(2 * y,     2 * y);
      K(index + 1, index)     += localK(2 * y + 1, 2 * y);
      K(index,     index + 1) += localK(2 * y,     2 * y + 1);
      K(index + 1, index + 1) += localK(2 * y + 1, 2 * y + 1);
    }

    // fill in the off-diagonal entries
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 4; i++)
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
void SQUARE_MESH::computeStiffnessMatrixSparse(SPARSE_MATRIX& K)
{
  TIMER functionTimer(__FUNCTION__);
  // allocate and clear
  if (K.rows() != _DOFs || K.cols() != _DOFs)
    K = SPARSE_MATRIX(_DOFs, _DOFs);
  K.setZero();

  typedef Eigen::Triplet<Real> TRIPLET;
  vector<TRIPLET> triplets;
  triplets.reserve(_DOFs * 4);

  for (int x = 0; x < _squares.size(); x++)
  {
    //MATRIX localK = _squares[x].computeForceJacobian();
    MATRIX localK = _squares[x].computeForceJacobianFast();

    // get the lookup indices
    int indices[4];
    for (int y = 0; y < 4; y++)
    {
      if (_vertexToIndex.find(_squares[x].vertex(y)) == _vertexToIndex.end())
        indices[y] = -1;
      else
        indices[y] = _vertexToIndex[_squares[x].vertex(y)];
    }

    // fill in the diagonal entries
    for (int y = 0; y < 4; y++)
    {
      if (indices[y] == -1) continue;

      int index = indices[y];
      /*
      K(index,     index)     += localK(2 * y,     2 * y);
      K(index + 1, index)     += localK(2 * y + 1, 2 * y);
      K(index,     index + 1) += localK(2 * y,     2 * y + 1);
      K(index + 1, index + 1) += localK(2 * y + 1, 2 * y + 1);
      */

      for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
          TRIPLET t(index + i, index + j, localK(2 * y + i, 2 * y + j));
          triplets.push_back(t);
        }
    }

    // fill in the off-diagonal entries
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 4; i++)
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
        /*
        K(rowStart,     colStart)     += localK(2 * i,     2 * j);
        K(rowStart + 1, colStart)     += localK(2 * i + 1, 2 * j);
        K(rowStart,     colStart + 1) += localK(2 * i,     2 * j + 1);
        K(rowStart + 1, colStart + 1) += localK(2 * i + 1, 2 * j + 1);
        */
      }
  }
  K.setFromTriplets(triplets.begin(), triplets.end());

  /*
  MATRIX fullK(K);
  EigenSolver<MATRIX> eigensolver(fullK);

  VECTOR realEigenvalues = -1 * eigensolver.eigenvalues().real();
  //cout << " K eigenvalues: " << eigensolver.eigenvalues().transpose() << endl;
  cout << " K eigenvalues: " << realEigenvalues.transpose() << endl;

  Real imaginary = eigensolver.eigenvalues().imag().norm();
  cout << " imaginary magnitude: " << imaginary << endl;

  int negatives = 0;
  for (unsigned int x = 0; x < realEigenvalues.size(); x++)
    if (realEigenvalues[x] < 0)
      negatives++;
  cout << " Negative eigenvalues: " << negatives << endl;
  */

  /*
  if (fabs(imaginary) > 1.0 || negatives > 0)
  {
    cout << " K: " << endl;
    cout << fullK << endl;
    exit(0);
  }
  */
}

///////////////////////////////////////////////////////////////////////
// compute the stiffness matrix
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::computeInvertibleStiffnessMatrixSparse(SPARSE_MATRIX& K)
{
  TIMER functionTimer(__FUNCTION__);
  // allocate and clear
  if (K.rows() != _DOFs || K.cols() != _DOFs)
    K = SPARSE_MATRIX(_DOFs, _DOFs);
  K.setZero();

  typedef Eigen::Triplet<Real> TRIPLET;
  vector<TRIPLET> triplets;
  triplets.reserve(_DOFs * 4);

  for (int x = 0; x < _squares.size(); x++)
  {
    //MATRIX localK = _squares[x].computeForceJacobianFast();
    MATRIX localK = _squares[x].computeInvertibleForceJacobian();

    // get the lookup indices
    int indices[4];
    for (int y = 0; y < 4; y++)
    {
      if (_vertexToIndex.find(_squares[x].vertex(y)) == _vertexToIndex.end())
        indices[y] = -1;
      else
        indices[y] = _vertexToIndex[_squares[x].vertex(y)];
    }

    // fill in the diagonal entries
    for (int y = 0; y < 4; y++)
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
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 4; i++)
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
// validate that the stiffness matrix does in fact yield the force
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::validateStiffness(const MATRIX& K, const VECTOR& f)
{
  VECTOR linearized = K * _u;

  VECTOR diff = f - linearized;

  cout << " f: " << endl;
  cout << f << endl;
  cout << " linear: " << endl;
  cout << linearized  << endl;
  cout << " u: " << endl;
  cout << _u << endl;
  cout << " square thinks displacement is: " << endl;
  cout << _squares[0].getDisplacement() << endl;
  cout << " diff: " << diff.dot(diff) << endl;

  //exit(0);
}

///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the stretch test
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::stepStretchForceTest(const Real stretch)
{
  for (int x = 0; x < _springConstraints.size(); x++)
    _springConstraints[x].first[0] += stretch;
  /*
  for (int x = 0; x < _constrainedVertices.size() / 2; x++)
  {
    //int left = _constrainedVertices[2 * x];
    int right = _constrainedVertices[2 * x + 1];

    //Real stretch = 0.01;
    //Real stretch = 0.001;

    //_vertices[left][0] -= stretch;
    _vertices[right][0] += stretch;
    //_vertices[right][0] -= stretch;
  }
  */
}

///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the stretch test
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::stepStretchTest(const Real stretch)
{
  for (int x = 0; x < _constrainedVertices.size(); x++)
  {
    int right = _constrainedVertices[x];
    if (_restVertices[right][0] > 0.5)
      _vertices[right][0] += stretch;
  }
}

///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the stretch test
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::stepCycleTest(const Real time)
{
  //Real amplitude = 0.1;
  Real amplitude = 0.05;
  Real amount = amplitude * cos(time - (M_PI / 8.0 + M_PI / 4.0) * 0.5);

  for (int x = 0; x < _constrainedVertices.size(); x++)
  {
    int right = _constrainedVertices[x];
    if (_restVertices[right][0] > 0.5)
      _vertices[right][0] += amount;
  }
}

///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the squash test
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::stepSquashForceTest(const Real stretch)
{
  for (int x = 0; x < _springConstraints.size(); x++)
    _springConstraints[x].first[0] -= stretch;
}

///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the squash test
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::stepSpinTest(const Real squash)
{
  static int frame = 0;
  for (int x = 0; x < _constrainedVertices.size(); x++)
  {
    int right = _constrainedVertices[x];
    if (_restVertices[right][0] > 0.1)
    {
      if (frame < 200)
        _vertices[right][0] -= squash;
      else if (frame < 400)
        _vertices[right][1] += squash;
      else if (frame < 750)
        _vertices[right][0] += squash;
      else if (frame < 1100)
        _vertices[right][1] -= squash;
      else if (frame < 1450)
        _vertices[right][0] -= squash;
      else if (frame < 1600)
        _vertices[right][1] += squash;
      else if (frame < 1800)
        _vertices[right][0] += squash;
    }
  }
  frame++;
  cout << " frame: " << frame << endl;
  /*
  for (int x = 0; x < _constrainedVertices.size() / 2; x++)
  {
    //int left = _constrainedVertices[2 * x];
    int right = _constrainedVertices[2 * x + 1];

    //Real stretch = 0.001;

    //_vertices[left][0] -= stretch;
    _vertices[right][0] -= squash;
    //_vertices[right][0] -= stretch;
  }
  */
}

///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the squash test
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::stepSquashTest(const Real squash)
{
  for (int x = 0; x < _constrainedVertices.size(); x++)
  {
    int right = _constrainedVertices[x];
    if (_restVertices[right][0] > 0.1)
      _vertices[right][0] -= squash;
  }
  /*
  for (int x = 0; x < _constrainedVertices.size() / 2; x++)
  {
    //int left = _constrainedVertices[2 * x];
    int right = _constrainedVertices[2 * x + 1];

    //Real stretch = 0.001;

    //_vertices[left][0] -= stretch;
    _vertices[right][0] -= squash;
    //_vertices[right][0] -= stretch;
  }
  */
}

///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the squash test
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::stepShearForceTest(const Real stretch)
{
  for (int x = 0; x < _springConstraints.size(); x++)
    _springConstraints[x].first[1] -= stretch;
}

///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the shear test
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::stepShearTest(const Real stretch)
{
  for (int x = 0; x < _constrainedVertices.size(); x++)
  {
    int right = _constrainedVertices[x];
    if (_restVertices[right][0] > 0.1)
      _vertices[right][1] -= stretch;
  }
  /*
  for (int x = 0; x < _constrainedVertices.size() / 2; x++)
  {
    //int left = _constrainedVertices[2 * x];
    int right = _constrainedVertices[2 * x + 1];

    //Real shear = 0.01;
    //Real shear = 0.1;

    //_vertices[left][0] -= stretch;
    _vertices[right][1] -= stretch;
    //_vertices[right][0] -= stretch;
  }
  */
}

///////////////////////////////////////////////////////////////////////
// get the average vertex position of the mesh
///////////////////////////////////////////////////////////////////////
VEC2 SQUARE_MESH::meshMean()
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
void SQUARE_MESH::boundingBox(VEC2& min, VEC2& max)
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
// get the area of the entire mesh
///////////////////////////////////////////////////////////////////////
Real SQUARE_MESH::meshArea() const
{
  Real area = 0;
  for (int x = 0; x < _squares.size(); x++)
    area += _squares[x].area();

  return area;
}

///////////////////////////////////////////////////////////////////////
// get the area of the smallest square
///////////////////////////////////////////////////////////////////////
Real SQUARE_MESH::smallestSquareArea() const
{
  Real area = fabs(_squares[0].area());
  for (int x = 0; x < _squares.size(); x++)
  {
    Real newest = _squares[x].area();
    if (fabs(newest) < area)
      area = fabs(newest);
  }

  return fabs(area);
}

///////////////////////////////////////////////////////////////////////
// get the area of the largest square
///////////////////////////////////////////////////////////////////////
Real SQUARE_MESH::largestSquareArea() const
{
  Real area = fabs(_squares[0].area());
  for (int x = 0; x < _squares.size(); x++)
  {
    Real newest = _squares[x].area();
    if (fabs(newest) > area)
      area = fabs(newest);
  }

  return fabs(area);
}

///////////////////////////////////////////////////////////////////////
// set a vertex position
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::setVertex(const VEC2& v, int index)
{
  assert(index < _vertices.size());

  _vertices[index] = v;
  uGather();
}

///////////////////////////////////////////////////////////////////////
// get the strain energy of the entire mesh
///////////////////////////////////////////////////////////////////////
Real SQUARE_MESH::psi() const
{
  Real final = 0;

  for (int x = 0; x < _squares.size(); x++)
    final += _squares[x].psi();

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute spring constraint forces
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::computeSpringForces()
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
void SQUARE_MESH::addSpringStiffnesses(SPARSE_MATRIX& K)
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
// probe the state (somehow)
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::probeState()
{
  cout << " Upper right psi: " << _squares[_upperRight].psi() << endl;
  cout << " Upper left psi: " << _squares[_upperLeft].psi() << endl;
  cout << " Lower right psi: " << _squares[_lowerRight].psi() << endl;
  
  //cout << " Upper right psi: " << _squares[_upperRight].computeForceVector() << endl;
  //cout << " Upper left psi: " << _squares[_upperLeft].computeForceVector() << endl;
  //cout << " Lower right psi: " << _squares[_lowerRight].computeForceVector() << endl;

  MATRIX ur = -1 * _squares[_upperRight].computeForceJacobianFast();
  MATRIX ul = -1 * _squares[_upperLeft].computeForceJacobianFast();
  MATRIX lr = -1 * _squares[_lowerRight].computeForceJacobianFast();

  // TK: BAD
  SelfAdjointEigenSolver<MATRIX> eigensolverur(ur);
  SelfAdjointEigenSolver<MATRIX> eigensolverul(ul);
  SelfAdjointEigenSolver<MATRIX> eigensolverlr(lr);
  cout << " upper right eigenvalues: " << eigensolverur.eigenvalues() << endl;
  cout << " upper left  eigenvalues: " << eigensolverul.eigenvalues() << endl;
  cout << " lower right eigenvalues: " << eigensolverlr.eigenvalues() << endl;
  //cout << " upper right Jacobian: " << _squares[_upperRight].computeForceJacobianFast() << endl;
  //cout << " upper left Jacobian: " << _squares[_upperLeft].computeForceJacobianFast() << endl;
  //cout << " lower right Jacobian: " << _squares[_lowerRight].computeForceJacobianFast() << endl;

  cout << " Upper right condition number: " << _squares[_upperRight].conditionNumber() << endl;
  cout << " Upper left condition number: " << _squares[_upperLeft].conditionNumber() << endl;
  cout << " Lower right condition number: " << _squares[_lowerRight].conditionNumber() << endl;
}

///////////////////////////////////////////////////////////////////////
// debugging function -- hand it the solution to 
// the compression test, see if it is recognized
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::setCompressionSolution(const Real squash)
{
  for (int x = 0; x < _unconstrainedVertices.size(); x++)
  {
    int index = _unconstrainedVertices[x];
    _vertices[index][0] = _restVertices[index][0] * squash;
  }
  uGather();
 
 /* 
  for (int x = 0; x < _constrainedVertices.size(); x++)
  {
    int index = _constrainedVertices[x];
    cout << " Constrained node " << x << ": " << _vertices[index].transpose() << endl;
  }
  */
}

///////////////////////////////////////////////////////////////////////
// copy the new inversion directions to the old ones
///////////////////////////////////////////////////////////////////////
void SQUARE_MESH::updateInversionDirections()
{
  for (unsigned int x = 0; x < _squares.size(); x++)
    _squares[x].updateInversionDirections();

  /*
  for (unsigned int x = 0; x < _squares.size(); x++)
    _squares[x].printReflectionInformation();
    */
}
