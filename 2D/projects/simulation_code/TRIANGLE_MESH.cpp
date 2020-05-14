#include "TRIANGLE_MESH.h"
#include "STVK.h"
#include <iostream>

#include <float.h>
#include <random>

using namespace std;

TRIANGLE_MESH::TRIANGLE_MESH(const Real poissonsRatio, const Real youngsModulus) : _DOFs(0)
{
  const Real E = youngsModulus;
  const Real nu = poissonsRatio;

  const Real lambda = (0.25)*( E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)));
  const Real mu = (0.25)*(E / (2.0 * (1 + nu)));
  _material = new STVK(lambda, mu);
}

TRIANGLE_MESH::~TRIANGLE_MESH()
{
  delete _material;
}

void TRIANGLE_MESH::buildBlob(const Real xPos, int sceneNum)
{
  _vertices.clear();
  _triangles.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();

  // build bottom (constrained) vertices
  for(int i = 0; i < 5; i++)
  {
    VEC2 v0(xPos + 0.1*i, -0.35);
    _vertices.push_back(v0);
    _restVertices.push_back(v0);
    _constrainedVertices.push_back(i);
  }

  // build tops to bottom triangles
  int vCount = 5;
  float scale = 5;
  float base = xPos + 0.05;
  for(int j = 4; j > 0; j--)
  {
    for(int i = 0; i < j; i++)
    {
      VEC2 v0(base + i*0.1, (sqrt(3.0) / 20.0)*(scale-j) -0.35 );
      _vertices.push_back(v0);
      _restVertices.push_back(v0);
      (j == 1)? _constrainedVertices.push_back(vCount + i) : _unconstrainedVertices.push_back(vCount + i);
    }
    vCount += j;
    base += 0.05;
  }

  int endofBaseVerts = vCount;

  // build left curved side vertices
  for(int i = 0; i < 3; i++)
  {
    float y;
    y = (i == 2)? (sqrt(3.0) / 20.0) + (sqrt(3.0) / 40.0) : (sqrt(3.0) / 40.0) + (i*(sqrt(3.0) / 20.0)/(i+1));
    VEC2 vo(xPos - 0.03 + (i+i)*0.01, y -0.35 );
    _vertices.push_back(vo);
    _restVertices.push_back(vo);
    _unconstrainedVertices.push_back(vCount);
    vCount++;
  }

  VEC2 vo(xPos + 0.05, (sqrt(3.0) / 20.0)*2 -0.35);
  _vertices.push_back(vo);
  _restVertices.push_back(vo);
  _unconstrainedVertices.push_back(vCount);
  vCount++;

  // build right curved side vertices
  for(int i = 0; i < 3; i++)
  {
    float y;
    y = (i == 2)? (sqrt(3.0) / 20.0) + (sqrt(3.0) / 40.0) : (sqrt(3.0) / 40.0) + (i*(sqrt(3.0) / 20.0)/(i+1));
    VEC2 vo2(xPos + 0.4 + 0.03 - (i+i)*0.01, y -0.35);
    _vertices.push_back(vo2);
    _restVertices.push_back(vo2);
    _unconstrainedVertices.push_back(vCount);
    vCount++;
  }

  VEC2 vo2(xPos + 0.35, (sqrt(3.0) / 20.0)*2 -0.35);
  _vertices.push_back(vo2);
  _restVertices.push_back(vo2);
  _unconstrainedVertices.push_back(vCount);
  vCount++;

  // build the base triangles
  int oldCount = 0;
  int newCount = 5;
  for(int max = 4; max > 0; max--)
  {
  for(int i = 0; i < max; i++)
  {
      vector<VEC2*> v(4);
      v[0] = &_vertices[oldCount + i];
      v[1] = &_vertices[oldCount + i + 1];
      v[2] = &_vertices[newCount + i];
      if(i < max - 1)
        v[3] = &_vertices[newCount + 1 + i];

      //bottom
      vector<VEC2*> triangle;
      triangle.push_back(v[0]);
      triangle.push_back(v[1]);
      triangle.push_back(v[2]);
      _triangles.push_back(TRIANGLE(_material, triangle));

      //top
      if(i < max - 1)
      {
        vector<VEC2*> top;
        top.push_back(v[2]);
        top.push_back(v[3]);
        top.push_back(v[1]);
        _triangles.push_back(TRIANGLE(_material, top));
      }
    }
    oldCount = newCount;
    newCount += max;
  }

  // build curved vertices
  for(int i = 0; i < 2; i++)
  {
    vector<VEC2*> v(12);
    v[0] = &_vertices[5 + i*3];
    v[1] = &_vertices[0 + i*4];
    v[2] = &_vertices[endofBaseVerts + i*4];
    v[3] = &_vertices[endofBaseVerts + 1 + i*4];
    v[4] = &_vertices[endofBaseVerts + 2 + i*4];
    v[5] = &_vertices[endofBaseVerts + 3 + i*4];
    v[6] = &_vertices[9 + i*2];
    v[7] = &_vertices[endofBaseVerts - (3 - i)];

    vector<VEC2*> triangle;
    triangle.push_back(v[0]);
    triangle.push_back(v[1]);
    triangle.push_back(v[2]);
    _triangles.push_back(TRIANGLE(_material, triangle));
    triangle.clear();

    triangle.push_back(v[0]);
    triangle.push_back(v[2]);
    triangle.push_back(v[3]);
    _triangles.push_back(TRIANGLE(_material, triangle));
    triangle.clear();

    triangle.push_back(v[0]);
    triangle.push_back(v[3]);
    triangle.push_back(v[4]);
    _triangles.push_back(TRIANGLE(_material, triangle));
    triangle.clear();

    triangle.push_back(v[0]);
    triangle.push_back(v[6]);
    triangle.push_back(v[4]);
    _triangles.push_back(TRIANGLE(_material, triangle));
    triangle.clear();

    triangle.push_back(v[4]);
    triangle.push_back(v[6]);
    triangle.push_back(v[5]);
    _triangles.push_back(TRIANGLE(_material, triangle));
    triangle.clear();

    triangle.push_back(v[5]);
    triangle.push_back(v[6]);
    triangle.push_back(v[7]);
    _triangles.push_back(TRIANGLE(_material, triangle));
  }

  // allocate the state vectors
  _DOFs = 2 * (_vertices.size() - _constrainedVertices.size());

  if (sceneNum)
    basisNoTranslation();
  else
  {
    setBasisReduction();
  }

  setMassMatrix();

  VECTOR zeros(_vertices.size()*2);
  VECTOR z2(_U.cols());

  z2.setZero();

  zeros.setZero();
  _u            = zeros;
  _q            = z2;
  _ra           = z2;
  _rv           = z2;
  _f            = z2;
  _fExternal    = zeros;
  _acceleration = zeros;
  _velocity     = zeros;

  // compute the reverse lookup
  computeVertexToIndexTable();
  createCoefs();
}

///////////////////////////////////////////////////////////////////////
// rebuild the vertex-to-index lookup
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeVertexToIndexTable()
{
  _vertexToIndex.clear();
  _allVertsToIndex.clear();
  for (unsigned int x = 0; x < _unconstrainedVertices.size(); x++)
  {
    VEC2* vertex = &_vertices[_unconstrainedVertices[x]];
    _vertexToIndex[vertex] = 2 * x;
    _allVertsToIndex[vertex] = 2 * x;
  }

  for (unsigned int x = 0; x < _constrainedVertices.size(); x++)
  {
    VEC2* vertex = &_vertices[_constrainedVertices[x]];
    _allVertsToIndex[vertex] = (2*_unconstrainedVertices.size()) + (2 * x);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::uScatter()
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
void TRIANGLE_MESH::uGather()
{
  for (unsigned int x = 0; x < _unconstrainedVertices.size(); x++)
  {
    int index = _unconstrainedVertices[x];
    _u[2 * x]     = _vertices[index][0] - _restVertices[index][0];
    _u[2 * x + 1] = _vertices[index][1] - _restVertices[index][1];
  }
}

void TRIANGLE_MESH::setMassMatrix()
{
  // matrix of size 2Nx2N
  MATRIX M(_vertices.size()*2,_vertices.size()*2);

  //for now we will make every vertex has mass 1
  M.setIdentity();
  M = M*15;
  M = M*_U;
  M = _U.transpose() * M;
  _mass = M;
}

void TRIANGLE_MESH::basisNoTranslation()
{
  MATRIX U(46,9);
  MATRIX svddiag(46, 9);

  U << 0.011446,  -0.002633,  -0.016704,   0.032098,  -0.009118,  -0.054865,  -0.054865,    0.088559,   0.067614,
     0.005271,  -0.001795,  -0.006527,   0.005844,  -0.015733,  -0.019406,  -0.019406,    0.006542,    -0.021152,
     0.004064,  -0.000470,  -0.017430,   0.022485,  -0.002303,  -0.085461,  -0.055144,    0.085461,     0.059172,
     0.009012,  -0.001448,  -0.000137,   0.003981,  -0.011465,  -0.000642,  -0.014926,    0.000642,    -0.017341,
    -0.004064,   0.000470,  -0.022485,   0.017430,   0.002303,  -0.085461,  -0.059079,    0.085461,     0.055184,
     0.009012,  -0.001448,   0.003981,  -0.000137,  -0.011465,   0.000642,  -0.017255,   -0.000642,    -0.014943,
    -0.011446,   0.002633,  -0.032098,   0.016704,   0.009118,  -0.088559,  -0.067503,    0.088559,     0.054882,
     0.005271,  -0.001795,   0.005844,  -0.006527,  -0.015733,   0.006542,  -0.020597,   -0.006542,    -0.019407,
     0.012217,  -0.003448,  -0.042263,   0.064138,  -0.005182,  -0.134988,  -0.079721,    0.134988,     0.101109,
     0.019461,  -0.005735,  -0.006665,   0.013395,  -0.022996,  -0.001420,  -0.013228,    0.001420,    -0.019055,
     0.000000,   0.000000,  -0.048845,   0.048845,   0.000000,  -0.138178,  -0.088842,    0.138178,     0.088817,
     0.022649,  -0.002830,   0.004947,   0.004947,  -0.019589,   0.000000,  -0.014737,   -0.000000,    -0.014741,
    -0.012217,   0.003448,  -0.064138,   0.042263,   0.005182,  -0.134988,  -0.101219,    0.134988,     0.079705,
     0.019461,  -0.005735,   0.013395,  -0.006665,  -0.022996,   0.001420,  -0.018808,   -0.001420,    -0.013217,
     0.012745,  -0.008094,  -0.085459,   0.113785,   0.002792,  -0.134270,  -0.063082,    0.134270,     0.086800,
     0.039573,  -0.011360,  -0.004205,   0.013985,  -0.020846,   0.013096,   0.000750,   -0.013096,    -0.009037,
    -0.012745,   0.008094,  -0.113785,   0.085459,  -0.002792,  -0.134270,  -0.086941,    0.134270,     0.063019,
     0.039573,  -0.011360,   0.013985,  -0.004205,  -0.020846,  -0.013096,  -0.008899,    0.013096,     0.000746,
     0.008654,  -0.002118,  -0.012365,   0.025727,  -0.019616,  -0.066554,  -0.030756,    0.066554,     0.069918,
     0.003497,  -0.001276,  -0.014071,   0.008843,  -0.032046,  -0.041551,  -0.032057,    0.041551,    -0.002250,
     0.011032,  -0.002679,  -0.019256,   0.036220,  -0.017551,  -0.098034,  -0.055586,    0.098034,     0.091227,
     0.004048,  -0.001544,  -0.015372,   0.010157,  -0.034182,  -0.037453,  -0.036060,    0.037453,    -0.005837,
     0.015469,  -0.003856,  -0.029304,   0.052425,  -0.017910,  -0.124694,  -0.070971,    0.124694,     0.114063,
     0.003343,  -0.001598,  -0.017875,   0.008649,  -0.034542,  -0.031420,  -0.033724,    0.031420,    -0.013374,
     0.014198,  -0.003726,  -0.039106,   0.064256,  -0.010398,  -0.140364,  -0.081399,    0.140364,     0.113799,
     0.010418,  -0.003566,  -0.013914,   0.010963,  -0.032107,  -0.019283,  -0.025475,    0.019283,    -0.013915,
    -0.008654,   0.002118,  -0.025727,   0.012365,   0.019616,  -0.066554,  -0.072897,    0.066554,     0.030760,
     0.003498,  -0.001276,   0.008843,  -0.014071,  -0.032046,   0.041551,  -0.005429,   -0.041551,    -0.032057,
    -0.011032,   0.002679,  -0.036220,   0.019256,   0.017551,  -0.098034,  -0.091957,    0.098034,     0.055597,
     0.004048,  -0.001544,   0.010157,  -0.015372,  -0.034182,   0.037453,  -0.006969,   -0.037453,    -0.036057,
    -0.015469,   0.003856,  -0.052425,   0.029304,   0.017910,  -0.124694,  -0.114055,    0.124694,     0.070974,
     0.003343,  -0.001598,   0.008649,  -0.017875,  -0.034542,   0.031420,  -0.012945,   -0.031420,    -0.033711,
    -0.014198,   0.003726,  -0.064256,   0.039106,   0.010398,  -0.140364,  -0.113929,    0.140364,     0.081385,
     0.010418,  -0.003566,   0.010963,  -0.013914,  -0.032107,   0.019283,  -0.013573,   -0.019283,    -0.025460,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000;

     JacobiSVD<MatrixXd> svd( U, ComputeFullV | ComputeFullU );
     svddiag = svd.matrixU();
     svddiag.conservativeResize(svddiag.rows(),9);
     _U = svddiag;

     // free space
     U.resize(0,0);
}

void TRIANGLE_MESH::setBasisReduction()
{
  MATRIX U(46,9);
  MATRIX T(46,2);
  MATRIX svddiag(46, 9);
  MATRIX intermediate(46, 11);
  svddiag.setIdentity();
  // matrix of deformations
  U << 0.011446,  -0.002633,  -0.016704,   0.032098,  -0.009118,  -0.054865,  -0.054865,    0.088559,   0.067614,
     0.005271,  -0.001795,  -0.006527,   0.005844,  -0.015733,  -0.019406,  -0.019406,    0.006542,    -0.021152,
     0.004064,  -0.000470,  -0.017430,   0.022485,  -0.002303,  -0.085461,  -0.055144,    0.085461,     0.059172,
     0.009012,  -0.001448,  -0.000137,   0.003981,  -0.011465,  -0.000642,  -0.014926,    0.000642,    -0.017341,
    -0.004064,   0.000470,  -0.022485,   0.017430,   0.002303,  -0.085461,  -0.059079,    0.085461,     0.055184,
     0.009012,  -0.001448,   0.003981,  -0.000137,  -0.011465,   0.000642,  -0.017255,   -0.000642,    -0.014943,
    -0.011446,   0.002633,  -0.032098,   0.016704,   0.009118,  -0.088559,  -0.067503,    0.088559,     0.054882,
     0.005271,  -0.001795,   0.005844,  -0.006527,  -0.015733,   0.006542,  -0.020597,   -0.006542,    -0.019407,
     0.012217,  -0.003448,  -0.042263,   0.064138,  -0.005182,  -0.134988,  -0.079721,    0.134988,     0.101109,
     0.019461,  -0.005735,  -0.006665,   0.013395,  -0.022996,  -0.001420,  -0.013228,    0.001420,    -0.019055,
     0.000000,   0.000000,  -0.048845,   0.048845,   0.000000,  -0.138178,  -0.088842,    0.138178,     0.088817,
     0.022649,  -0.002830,   0.004947,   0.004947,  -0.019589,   0.000000,  -0.014737,   -0.000000,    -0.014741,
    -0.012217,   0.003448,  -0.064138,   0.042263,   0.005182,  -0.134988,  -0.101219,    0.134988,     0.079705,
     0.019461,  -0.005735,   0.013395,  -0.006665,  -0.022996,   0.001420,  -0.018808,   -0.001420,    -0.013217,
     0.012745,  -0.008094,  -0.085459,   0.113785,   0.002792,  -0.134270,  -0.063082,    0.134270,     0.086800,
     0.039573,  -0.011360,  -0.004205,   0.013985,  -0.020846,   0.013096,   0.000750,   -0.013096,    -0.009037,
    -0.012745,   0.008094,  -0.113785,   0.085459,  -0.002792,  -0.134270,  -0.086941,    0.134270,     0.063019,
     0.039573,  -0.011360,   0.013985,  -0.004205,  -0.020846,  -0.013096,  -0.008899,    0.013096,     0.000746,
     0.008654,  -0.002118,  -0.012365,   0.025727,  -0.019616,  -0.066554,  -0.030756,    0.066554,     0.069918,
     0.003497,  -0.001276,  -0.014071,   0.008843,  -0.032046,  -0.041551,  -0.032057,    0.041551,    -0.002250,
     0.011032,  -0.002679,  -0.019256,   0.036220,  -0.017551,  -0.098034,  -0.055586,    0.098034,     0.091227,
     0.004048,  -0.001544,  -0.015372,   0.010157,  -0.034182,  -0.037453,  -0.036060,    0.037453,    -0.005837,
     0.015469,  -0.003856,  -0.029304,   0.052425,  -0.017910,  -0.124694,  -0.070971,    0.124694,     0.114063,
     0.003343,  -0.001598,  -0.017875,   0.008649,  -0.034542,  -0.031420,  -0.033724,    0.031420,    -0.013374,
     0.014198,  -0.003726,  -0.039106,   0.064256,  -0.010398,  -0.140364,  -0.081399,    0.140364,     0.113799,
     0.010418,  -0.003566,  -0.013914,   0.010963,  -0.032107,  -0.019283,  -0.025475,    0.019283,    -0.013915,
    -0.008654,   0.002118,  -0.025727,   0.012365,   0.019616,  -0.066554,  -0.072897,    0.066554,     0.030760,
     0.003498,  -0.001276,   0.008843,  -0.014071,  -0.032046,   0.041551,  -0.005429,   -0.041551,    -0.032057,
    -0.011032,   0.002679,  -0.036220,   0.019256,   0.017551,  -0.098034,  -0.091957,    0.098034,     0.055597,
     0.004048,  -0.001544,   0.010157,  -0.015372,  -0.034182,   0.037453,  -0.006969,   -0.037453,    -0.036057,
    -0.015469,   0.003856,  -0.052425,   0.029304,   0.017910,  -0.124694,  -0.114055,    0.124694,     0.070974,
     0.003343,  -0.001598,   0.008649,  -0.017875,  -0.034542,   0.031420,  -0.012945,   -0.031420,    -0.033711,
    -0.014198,   0.003726,  -0.064256,   0.039106,   0.010398,  -0.140364,  -0.113929,    0.140364,     0.081385,
     0.010418,  -0.003566,   0.010963,  -0.013914,  -0.032107,   0.019283,  -0.013573,   -0.019283,    -0.025460,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,    0.000000,     0.000000;

     T << 1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000,
          1.000000,  0.000000,
          0.000000,  1.000000;

  T = (1.0/pow(23.0, 0.5))*T;

  JacobiSVD<MatrixXd> svd( U, ComputeFullV | ComputeFullU );
  svddiag = svd.matrixU();
  svddiag.conservativeResize(svddiag.rows(),9);

  intermediate.col(0) = T.col(0);
  intermediate.col(1) = T.col(1);
  for(int i = 2; i < 11; i++)
  {
    intermediate.col(i) = svddiag.col(i - 2) - (T.col(0).transpose()*svddiag.col(i - 2))*T.col(0);
    intermediate.col(i) = intermediate.col(i) - (T.col(1).transpose()*intermediate.col(i))*T.col(1);
  }

  _U = intermediate;

  U.resize(0,0);
  T.resize(0,0);
}

void TRIANGLE_MESH::qTou()
{
  _u = _U * _q;
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

void TRIANGLE_MESH::addSingleForce(const VEC2& bodyForce, int vertex)
{
  std::vector<int>::iterator it = std::find(_constrainedVertices.begin(), _constrainedVertices.end(), vertex);
  int velocity_index;
  if(it != _constrainedVertices.end())
  {
    velocity_index = std::distance(_constrainedVertices.begin(), it);
    _fExternal[_unconstrainedVertices.size()*2 + velocity_index*2]     += bodyForce[0];
    _fExternal[_unconstrainedVertices.size()*2 + velocity_index*2 + 1] += bodyForce[1];
  }
  else
  {
    it = std::find(_unconstrainedVertices.begin(), _unconstrainedVertices.end(), vertex);
    velocity_index = std::distance(_unconstrainedVertices.begin(), it);
    _fExternal[velocity_index*2]     += bodyForce[0];
    _fExternal[velocity_index*2 + 1] += bodyForce[1];
  }
}
///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the shear test
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::stepShearTest(const Real shear)
{
  for (unsigned int x = 0; x < _constrainedVertices.size(); x++)
  {
    int right = _constrainedVertices[x];
    if (_restVertices[right][1] > -.35)
      _vertices[right][0] += shear;
  }
}

///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the stretch test
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::stepStretchTest(const Real stretch)
{
  for (unsigned int x = 0; x < _constrainedVertices.size(); x++)
  {
    int right = _constrainedVertices[x];
    if (_restVertices[right][0] > 0.15)
      _vertices[right][0] += stretch;
  }
}

void TRIANGLE_MESH::stretch2(const Real stretch)
{
  for (unsigned int x = 0; x < _constrainedVertices.size(); x++)
  {
    int right = _constrainedVertices[x];
    if (_restVertices[right][1] > -.35)
      _vertices[right][1] += stretch;
  }
}


void TRIANGLE_MESH::createCoefs()
{
  MATRIX m(_vertices.size()*2,_vertices.size()*2);
  MATRIX constTemp(_vertices.size()*2,_vertices.size()*2);
  TENSOR3 cubq(_vertices.size()*2,_vertices.size()*2, _vertices.size()*2);
  TENSOR3 quadl(_vertices.size()*2,_vertices.size()*2, _vertices.size()*2);
  TENSOR4 tempQuad(_vertices.size()*2,_vertices.size()*2, _vertices.size()*2,_vertices.size()*2);
  TENSOR4 temp(_vertices.size()*2,_vertices.size()*2, _vertices.size()*2,_vertices.size()*2);
  m.setZero();
  constTemp.setZero();

  int n_of_triangles = _triangles.size();

  for(int x = 0; x < n_of_triangles; x++)
  {
    //get our current triangle
    TRIANGLE current = getTriangle(x);

    //get the const coef
    MATRIX constCoef = current.getConst();
    MATRIX linearCoef = current.getLinear();
    TENSOR3 cubquad = current.getQuadCubic();
    TENSOR3 quadlin = current.getLinearQuad();
    TENSOR4 quadCoef = current.getQuad();
    TENSOR4 cubicCoef = current.getCubic();

    //add the forces into the right place in the global linear coef
    for(int i = 0; i < 3; i++)
    {
      //get current triangle vertex y
      VEC2* current_iv = _triangles[x].vertex(i);

      //if the vertex we're on is unconstrained, add it to the global linear coef
      if(_vertexToIndex.find(current_iv) != _vertexToIndex.end())
      {
        //find the global index using the vertexToIndex map
        int global_index_i = _vertexToIndex[current_iv];

        for(int j = 0; j < 3; j++)
        {
          bool stiffness = false;
          VEC2* current_jv = _triangles[x].vertex(j);

          // set values in linear coefficient of cubic polynomial
          int global_index_j = _allVertsToIndex[current_jv];
          m(global_index_i, global_index_j) += linearCoef(2*i, 2*j);
          m(global_index_i, global_index_j + 1) += linearCoef(2*i, 2*j + 1);
          m(global_index_i + 1, global_index_j) += linearCoef(2*i + 1, 2*j);
          m(global_index_i + 1, global_index_j + 1) += linearCoef(2*i + 1, 2*j + 1);

          if(_vertexToIndex.find(current_jv) != _vertexToIndex.end())
          {
            stiffness = true;

            // set values of constant coefficient in quadratic polynomial
            constTemp(global_index_i, global_index_j) += constCoef(2*i, 2*j);
            constTemp(global_index_i, global_index_j + 1) += constCoef(2*i, 2*j + 1);
            constTemp(global_index_i + 1, global_index_j) += constCoef(2*i + 1, 2*j);
            constTemp(global_index_i + 1, global_index_j + 1) += constCoef(2*i + 1, 2*j + 1);
          }

          for(int k = 0; k < 3; k++)
          {
            VEC2* current_kv = _triangles[x].vertex(k);

            int global_index_k = _allVertsToIndex[current_kv];

            // set values in quadratic coefficient in cubic polynomial
            cubq._tensor[global_index_k](global_index_i, global_index_j) += cubquad._tensor[2*k](2*i, 2*j);
            cubq._tensor[global_index_k](global_index_i, global_index_j + 1) += cubquad._tensor[2*k](2*i, 2*j + 1);

            cubq._tensor[global_index_k](global_index_i + 1, global_index_j) += cubquad._tensor[2*k](2*i + 1, 2*j);
            cubq._tensor[global_index_k](global_index_i + 1, global_index_j + 1) += cubquad._tensor[2*k](2*i + 1, 2*j + 1);

            cubq._tensor[global_index_k + 1](global_index_i, global_index_j) += cubquad._tensor[2*k + 1](2*i, 2*j);
            cubq._tensor[global_index_k + 1](global_index_i, global_index_j + 1) += cubquad._tensor[2*k + 1](2*i, 2*j + 1);

            cubq._tensor[global_index_k + 1](global_index_i + 1, global_index_j) += cubquad._tensor[2*k + 1](2*i + 1, 2*j);
            cubq._tensor[global_index_k + 1](global_index_i + 1, global_index_j + 1) += cubquad._tensor[2*k + 1](2*i + 1, 2*j + 1);

            if(stiffness)
            {
              // set values in linear coefficient in quadratic polynomial
              quadl._tensor[global_index_k](global_index_i, global_index_j) += quadlin._tensor[2*k](2*i, 2*j);
              quadl._tensor[global_index_k](global_index_i, global_index_j + 1) += quadlin._tensor[2*k](2*i, 2*j + 1);

              quadl._tensor[global_index_k](global_index_i + 1, global_index_j) += quadlin._tensor[2*k](2*i + 1, 2*j);
              quadl._tensor[global_index_k](global_index_i + 1, global_index_j + 1) += quadlin._tensor[2*k](2*i + 1, 2*j + 1);

              quadl._tensor[global_index_k + 1](global_index_i, global_index_j) += quadlin._tensor[2*k + 1](2*i, 2*j);
              quadl._tensor[global_index_k + 1](global_index_i, global_index_j + 1) += quadlin._tensor[2*k + 1](2*i, 2*j + 1);

              quadl._tensor[global_index_k + 1](global_index_i + 1, global_index_j) += quadlin._tensor[2*k + 1](2*i + 1, 2*j);
              quadl._tensor[global_index_k + 1](global_index_i + 1, global_index_j + 1) += quadlin._tensor[2*k + 1](2*i + 1, 2*j + 1);
            }

            for(int l = 0; l < 3; l++)
            {
              VEC2* current_lv = _triangles[x].vertex(l);

              int global_index_l = _allVertsToIndex[current_lv];

              // set coefficients in cubic term of cubic polynomial
              temp._tensor[global_index_l]._tensor[global_index_k](global_index_i, global_index_j) += cubicCoef._tensor[2*l]._tensor[2*k](2*i, 2*j);
              temp._tensor[global_index_l + 1]._tensor[global_index_k](global_index_i, global_index_j) += cubicCoef._tensor[2*l + 1]._tensor[2*k](2*i, 2*j);

              temp._tensor[global_index_l]._tensor[global_index_k + 1](global_index_i, global_index_j) += cubicCoef._tensor[2*l]._tensor[2*k + 1](2*i, 2*j);
              temp._tensor[global_index_l + 1]._tensor[global_index_k + 1](global_index_i, global_index_j) += cubicCoef._tensor[2*l + 1]._tensor[2*k + 1](2*i, 2*j);

              temp._tensor[global_index_l]._tensor[global_index_k](global_index_i, global_index_j + 1) += cubicCoef._tensor[2*l]._tensor[2*k](2*i, 2*j + 1);
              temp._tensor[global_index_l + 1]._tensor[global_index_k](global_index_i, global_index_j + 1) += cubicCoef._tensor[2*l + 1]._tensor[2*k](2*i, 2*j + 1);
              temp._tensor[global_index_l]._tensor[global_index_k + 1](global_index_i, global_index_j + 1) += cubicCoef._tensor[2*l]._tensor[2*k + 1](2*i, 2*j + 1);
              temp._tensor[global_index_l + 1]._tensor[global_index_k + 1](global_index_i, global_index_j + 1) += cubicCoef._tensor[2*l + 1]._tensor[2*k + 1](2*i, 2*j + 1);

              temp._tensor[global_index_l]._tensor[global_index_k](global_index_i + 1, global_index_j) += cubicCoef._tensor[2*l]._tensor[2*k](2*i + 1, 2*j);
              temp._tensor[global_index_l + 1]._tensor[global_index_k](global_index_i + 1, global_index_j) += cubicCoef._tensor[2*l + 1]._tensor[2*k](2*i + 1, 2*j);

              temp._tensor[global_index_l]._tensor[global_index_k + 1](global_index_i + 1, global_index_j) += cubicCoef._tensor[2*l]._tensor[2*k + 1](2*i + 1, 2*j);
              temp._tensor[global_index_l + 1]._tensor[global_index_k + 1](global_index_i + 1, global_index_j) += cubicCoef._tensor[2*l + 1]._tensor[2*k + 1](2*i + 1, 2*j);

              temp._tensor[global_index_l]._tensor[global_index_k](global_index_i + 1, global_index_j + 1) += cubicCoef._tensor[2*l]._tensor[2*k](2*i + 1, 2*j + 1);
              temp._tensor[global_index_l + 1]._tensor[global_index_k](global_index_i + 1, global_index_j + 1) += cubicCoef._tensor[2*l + 1]._tensor[2*k](2*i + 1, 2*j + 1);
              temp._tensor[global_index_l]._tensor[global_index_k + 1](global_index_i + 1, global_index_j + 1) += cubicCoef._tensor[2*l]._tensor[2*k + 1](2*i + 1, 2*j + 1);
              temp._tensor[global_index_l + 1]._tensor[global_index_k + 1](global_index_i + 1, global_index_j + 1) += cubicCoef._tensor[2*l + 1]._tensor[2*k + 1](2*i + 1, 2*j + 1);

              if(stiffness)
              {
                // set coefficients in quadratic term of quadratic polynomial
                tempQuad._tensor[global_index_l]._tensor[global_index_k](global_index_i, global_index_j) += quadCoef._tensor[2*l]._tensor[2*k](2*i, 2*j);
                tempQuad._tensor[global_index_l + 1]._tensor[global_index_k](global_index_i, global_index_j) += quadCoef._tensor[2*l + 1]._tensor[2*k](2*i, 2*j);

                tempQuad._tensor[global_index_l]._tensor[global_index_k + 1](global_index_i, global_index_j) += quadCoef._tensor[2*l]._tensor[2*k + 1](2*i, 2*j);
                tempQuad._tensor[global_index_l + 1]._tensor[global_index_k + 1](global_index_i, global_index_j) += quadCoef._tensor[2*l + 1]._tensor[2*k + 1](2*i, 2*j);

                tempQuad._tensor[global_index_l]._tensor[global_index_k](global_index_i, global_index_j + 1) += quadCoef._tensor[2*l]._tensor[2*k](2*i, 2*j + 1);
                tempQuad._tensor[global_index_l + 1]._tensor[global_index_k](global_index_i, global_index_j + 1) += quadCoef._tensor[2*l + 1]._tensor[2*k](2*i, 2*j + 1);
                tempQuad._tensor[global_index_l]._tensor[global_index_k + 1](global_index_i, global_index_j + 1) += quadCoef._tensor[2*l]._tensor[2*k + 1](2*i, 2*j + 1);
                tempQuad._tensor[global_index_l + 1]._tensor[global_index_k + 1](global_index_i, global_index_j + 1) += quadCoef._tensor[2*l + 1]._tensor[2*k + 1](2*i, 2*j + 1);

                tempQuad._tensor[global_index_l]._tensor[global_index_k](global_index_i + 1, global_index_j) += quadCoef._tensor[2*l]._tensor[2*k](2*i + 1, 2*j);
                tempQuad._tensor[global_index_l + 1]._tensor[global_index_k](global_index_i + 1, global_index_j) += quadCoef._tensor[2*l + 1]._tensor[2*k](2*i + 1, 2*j);

                tempQuad._tensor[global_index_l]._tensor[global_index_k + 1](global_index_i + 1, global_index_j) += quadCoef._tensor[2*l]._tensor[2*k + 1](2*i + 1, 2*j);
                tempQuad._tensor[global_index_l + 1]._tensor[global_index_k + 1](global_index_i + 1, global_index_j) += quadCoef._tensor[2*l + 1]._tensor[2*k + 1](2*i + 1, 2*j);

                tempQuad._tensor[global_index_l]._tensor[global_index_k](global_index_i + 1, global_index_j + 1) += quadCoef._tensor[2*l]._tensor[2*k](2*i + 1, 2*j + 1);
                tempQuad._tensor[global_index_l + 1]._tensor[global_index_k](global_index_i + 1, global_index_j + 1) += quadCoef._tensor[2*l + 1]._tensor[2*k](2*i + 1, 2*j + 1);
                tempQuad._tensor[global_index_l]._tensor[global_index_k + 1](global_index_i + 1, global_index_j + 1) += quadCoef._tensor[2*l]._tensor[2*k + 1](2*i + 1, 2*j + 1);
                tempQuad._tensor[global_index_l + 1]._tensor[global_index_k + 1](global_index_i + 1, global_index_j + 1) += quadCoef._tensor[2*l + 1]._tensor[2*k + 1](2*i + 1, 2*j + 1);
              }

            }
          }
        }
      }
    }

    constCoef.resize(0,0);
    linearCoef.resize(0,0);
    quadCoef.clear();
    quadlin.clear();
    cubicCoef.clear();
    cubquad.clear();
  }

  MATRIX transpose = _U.transpose();

  // set reduced coefficients for the cubic polynomial

  // linear term
  _linearCoef = transpose*(m*_U);
  m.resize(0,0);

  // cubic coefficient
  _cubicCoef = temp.modeFourProduct(transpose);
  temp.clear();
  _cubicCoef = _cubicCoef.modeThreeProduct(transpose);
  _cubicCoef = _cubicCoef.modeTwoProduct(transpose);
  _cubicCoef = _cubicCoef.modeOneProduct(transpose);

  // quadratic coefficient of cubic polynomial
  _cubicquad = cubq.modeThreeProduct(transpose);
  cubq.clear();
  _cubicquad = _cubicquad.modeTwoProduct(transpose);
  _cubicquad = _cubicquad.modeOneProduct(transpose);

  // set reduced coefficients for quadratic polynomial

  // constant term
  _constCoef = transpose*(constTemp*_U);

  // linear term
  _quadlinear = quadl.modeThreeProduct(transpose);
  quadl.clear();
  _quadlinear = _quadlinear.modeTwoProduct(transpose);
  _quadlinear = _quadlinear.modeOneProduct(transpose);

  // quadratic term
  _quadraticCoef = tempQuad.modeFourProduct(transpose);
  tempQuad.clear();
  _quadraticCoef = _quadraticCoef.modeThreeProduct(transpose);
  _quadraticCoef = _quadraticCoef.modeTwoProduct(transpose);
  _quadraticCoef = _quadraticCoef.modeOneProduct(transpose);

  // free memory
  transpose.resize(0,0);
  temp.clear();
  cubq.clear();
  m.resize(0,0);

}
///////////////////////////////////////////////////////////////////////
//this will find the global force vector of forces on each unrestrained
//vertex.
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeMaterialForces()
{
  //the global vector will be v= [f_0, f_1 ...] where each f_i = [x,y] (column vectors)
  //and forces are only for unconstrained vertices. This means the size of this vector is

  VECTOR global_vector(_q.size());
  global_vector.setZero();

  //add linear term to global vector
  global_vector += (_linearCoef * _q);

  // compute the cubic term
  TENSOR3 temp = _cubicCoef.modeFourProduct(_q);
  MATRIX tempMatrix = temp.modeThreeProduct(_q);

  //add the cubic term to the global vector
  global_vector += (tempMatrix*_q);

  //compute the quadratic term, add to global vector
  tempMatrix = _cubicquad.modeThreeProduct(_q);
  global_vector += (tempMatrix * _q);

  // clear the memory of temporary things
  tempMatrix.resize(0,0);
  temp.clear();

  //store global vector into _f
  _f = global_vector;

}

void TRIANGLE_MESH::wackyCollision()
{
  float kw = 100; // spring constant of wall
  float l = 0.1; // dampening force constant

  for(unsigned int y = 0; y < _walls.size(); y++)
  {
    for(int x = 0; x < _vertices.size(); x++ )
    {
      float diffx;
      if(_walls[y].point()[0] > 0)
        diffx = _walls[y].point()[0] - (_vertices[x][0]);
      else
        diffx = _walls[y].point()[0] - (_vertices[x][0]);

      float diffy =  _walls[y].point()[1] - (_vertices[x][1]);
      int velocity_index;
      std::vector<int>::iterator it = std::find(_constrainedVertices.begin(), _constrainedVertices.end(), x);
      if(it != _constrainedVertices.end())
      {
        velocity_index = std::distance(_constrainedVertices.begin(), it);
        velocity_index =  _unconstrainedVertices.size()*2 + velocity_index*2;
      }
      else
      {
        it = std::find(_unconstrainedVertices.begin(), _unconstrainedVertices.end(), x);
        velocity_index = std::distance(_unconstrainedVertices.begin(), it);
        velocity_index = velocity_index*2;
      }

      if((diffx >= 0 && _walls[y].point()[0] < 0) || (diffx <= 0 && _walls[y].point()[0] > 0) ) //did it hit a side wall?
      {
        VEC2 force;
        force = ( kw * abs(diffx) * _walls[y].normal() ); //apply spring force of wall
        force += ( l * _velocity[velocity_index] * _walls[y].normal() ); //apply dampening force
        _fExternal[velocity_index] += force[0];
        _fExternal[velocity_index + 1] += force[1];
      }
      if(diffy >= 0 && _walls[y].point()[1] != 0) // did it hit the floor?
      {
        VEC2 force;
        force = (kw * diffy * _walls[y].normal() ); //apply spring force of wall
        force += ( l * _velocity[velocity_index + 1] * _walls[y].normal() ); //apply dampening force
        _fExternal[velocity_index] += force[0];
        _fExternal[velocity_index + 1] += force[1];
      }
    }
  }
}
// collision detection
void TRIANGLE_MESH::checkCollision()
{
  float kw = 100; // spring constant of wall
  float l = 0.1; // dampening force constant

  for(unsigned int y = 0; y < _walls.size(); y++)
  {
    for(int x = 0; x < _vertices.size(); x++ )
    {
      float diffx;
      if(_walls[y].point()[0] > 0)
        diffx = _walls[y].point()[0] - (_vertices[x][0]);
      else
        diffx = _walls[y].point()[0] - (_vertices[x][0]);

      float diffy =  _walls[y].point()[1] - (_vertices[x][1]);

      if((diffx >= 0 && _walls[y].point()[0] < 0) || (diffx <= 0 && _walls[y].point()[0] > 0) ) //did it hit a side wall?
      {
        addBodyForce( kw * abs(diffx) * _walls[y].normal() ); //apply spring force of wall
        addBodyForce( l * _velocity[2*x] * _walls[y].normal() ); //apply dampening force
        addSingleForce(abs(diffx) * _walls[y].normal(), x);
        break;
      }
      if(diffy >= 0 && _walls[y].point()[1] != 0) // did it hit the floor?
      {
          addBodyForce( kw * diffy * _walls[y].normal() ); //apply spring force of wall
          addBodyForce( l * _velocity[2*x + 1] * _walls[y].normal() ); //apply dampening force
          addSingleForce(abs(diffy) * _walls[y].normal(), x);
          break;
      }
    }
  }
}

// motion step using Euler Lagrange
void TRIANGLE_MESH::stepMotion(float dt, const VEC2& outerForce)
{
  //make stiffness Matrix K. size is 2*unrestrained vertices x  2*unrestrained vertices
  MATRIX K(_q.size(),_q.size() );
  MATRIX D(_q.size(),_q.size() );
  MATRIX inverse;
  float alpha = 0.01; // constant for damping
  float beta = 0.02;  // constant for damping

  checkCollision();

  // Newton Raphson Iteration, but j-max is 1 so no need to write the loop
  //step 1: compute K
  K.setZero();
  computeStiffnessMatrix(K);

  // step 2: compute D
  D = alpha*_mass - beta*K;


  // step 3: calculate f_external
  VECTOR reducedF = _U.transpose() * _fExternal;

  // step 4: compute R(q+1)
  computeMaterialForces();

  // step 5: calculate a1 - a6 with beta  0.25 and gamma = 0.5
  float betat = 0.25;
  float gamma = 0.5;
  float a1 = 1.0 / (betat* pow(dt, 2));
  float a2 = 1.0 / (betat * dt);
  float a3 = (1.0 - 2*betat) / (2*betat);
  float a4 = gamma / (betat*dt);
  float a5 = 1.0 - (gamma/betat);
  float a6 = (1.0 - (gamma/(2*betat)))*dt;

  // step 6: solve the equations
  VECTOR rightSolve = -1*((-1*a3*_mass + a6*D)*_ra + (-1*a2*_mass + a5*D)*_rv - _f - reducedF);
  MATRIX leftMatrix = a1*_mass + a4*D - K;
  inverse = leftMatrix.inverse().eval();

  VECTOR dq = inverse*rightSolve;

  // free space
  leftMatrix.resize(0,0);
  inverse.resize(0,0);

  // step 7: update q and u
  _q += dq;
  qTou();

  // step 8: update all node positions w new displacement vector
  int unconstrained = _unconstrainedVertices.size();
  uScatter();

  for(int x = 0; x < _constrainedVertices.size(); x++)
  {
    VEC2 displacement;
    displacement[0] = _u[2*unconstrained + 2*x];
    displacement[1] = _u[2*unconstrained + 2*x + 1];
    _vertices[_constrainedVertices[x]] = _restVertices[_constrainedVertices[x]]+ displacement;
  }

  // step 9: calculate velocity and accleration
  VECTOR newVel = a4*dq + a5*_rv + a6*_ra;
  _ra = a1*dq - a2*_rv - a3*_ra;
  _rv = newVel;
  _velocity = _U*_rv;

  _fExternal.setZero();
  _f.setZero();

}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::stepQuasistatic()
{
  //make stiffness Matrix K. size is 2*unrestrained vertices x  2*unrestrained vertices
  MATRIX K(_q.size(),_q.size());

  //step 1: compute K
  K.setZero();
  computeStiffnessMatrix(K);

  //step 2: compute internal material forces, R(uq)
  computeMaterialForces();

  // //step 3: external forces transform
  VECTOR reducedF = _U.transpose() * _fExternal;

  // step 4: form the residual (r = F + E)
  VECTOR r2 = -1*(_f + reducedF);

  //step 5: compute x = K .inverse().eval()  * r
  MATRIX inverse2 = K.inverse().eval();
  VECTOR x2 = inverse2*r2;

  //step 6: add solution x to displamcement vector _u
  _q += x2;
  qTou();

  //step 7: update all node positions w new displacement vector
  uScatter();

  //reset forces to 0
  _fExternal.setZero();
  _f.setZero();

  static int counter = 0;
  cout << " Quasistatic step: " << counter << endl;
  counter++;
  return true;
}

///////////////////////////////////////////////////////////////////////
// compute the stiffness matrix
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeStiffnessMatrix(MATRIX& K)
{
  //can assume K is the correct size, 2V x 2V

  // compute quadratic term and add to K
  TENSOR3 temp = _quadraticCoef.modeFourProduct(_q);
  K = temp.modeThreeProduct(_q);

  // compute linear term and add to K
  K += _quadlinear.modeThreeProduct(_q);

  // add constant term to K
  K += _constCoef;

  // free space
  temp.clear();
}
