#include "TRIANGLE_MESH.h"
// #include "STVK.h"
#include "NEOHOOKEAN.h"
#include <iostream>
#include <fstream>

#include <float.h>
#include <random>

using namespace std;

TRIANGLE_MESH::TRIANGLE_MESH(const Real poissonsRatio, const Real youngsModulus) : _DOFs(0)
{
  const Real E = youngsModulus;
  const Real nu = poissonsRatio;

  const Real lambda = (0.25)*( E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)));
  const Real mu = (0.25)*(E / (2.0 * (1 + nu)));
  // _material = new STVK(lambda, mu);
  _material = new NEOHOOKEAN(lambda, mu);
}

TRIANGLE_MESH::~TRIANGLE_MESH()
{
  delete _material;
}

void TRIANGLE_MESH::buildBlob(int sceneNum, const char* filename, bool create_basis, int cols)
{
  _vertices.clear();
  _triangles.clear();
  _constrainedVertices.clear();
  _unconstrainedVertices.clear();
  int vertCount;
  int size;

  // first we will get all the vertices and their coordinates. open the node file.
  string nodeFile = filename + string(".node");
  std::ifstream nodes(nodeFile);

  if (nodes.is_open())
  {
    std::string line;

    //get first line
    std::getline(nodes, line);

    // split first line
    std::istringstream iss(line);

    if(!(iss >> vertCount >> size)) // this line inside the if statement stores the values
    {
      cout << "error: issue with file format" << endl;
      nodes.close();
      exit(0);
    }

    // we are only dealing with 2d now, but maybe we want to keep the option open for 3d implementation
    if(size == 2)
    {
      VEC2 mins(1000,1000);
      VEC2 maxs(-1000,-1000);
      // store every vertex into our array
      for(int i = 0; i < vertCount; i++)
      {
        //get the next line
        std::getline(nodes, line);
        std::istringstream iss2(line);

        // x y coordinates to store
        int index;
        double x;
        double y;

        // if it doesnt match this format, there's an error
        if(!(iss2 >> index >> x >> y))
        {
          cout << "error: issue with file format" << endl;
          nodes.close();
          exit(0);
        }

        // store this new vertex
        VEC2 vert(x, y);
        _vertices.push_back(vert);
        _restVertices.push_back(vert);

        // find the max and min x,y positions
        if (vert[0] < mins[0])
          mins[0] = vert[0];
        if (vert[1] < mins[1])
          mins[1] = vert[1];
        if (vert[0] > maxs[0])
          maxs[0] = vert[0];
        if (vert[1] > maxs[1])
          maxs[1] = vert[1];

      }

      // shift vertices down so that the object rests on the floor, then constrain some vertices.
      for(int i = 0; i < vertCount; i++)
      {
        _vertices[i][1] += (-0.95 - mins[1]);
        _restVertices[i][1] += (-0.95 - mins[1]);

        (_restVertices[i][1] < -0.85 || (_restVertices[i][1] > maxs[1] - 1.05 - mins[1]))? _constrainedVertices.push_back(i) : _unconstrainedVertices.push_back(i);
      }
    }
  }
  else // error! shut down
  {
    cout << " Could not open file " << nodeFile << "!!!" << endl;
    exit(0);
  }

  // we are done with the nodes file! close it.
  nodes.close();

  // open the file with the triangle definitions and create our array of triangles
  nodeFile = filename + string(".ele");
  std::ifstream polys(nodeFile);

  if (polys.is_open())
  {
    std::string line;
    int numTri; // the number of triangles we are creating
    int numVert; // we need to make sure we are actually given 3 verts for a triangle

    //get first line
    std::getline(polys, line);

    // split first line
    std::istringstream iss(line);

    if(!(iss >> numTri >> numVert)) // this line inside the if statement stores the values
    {
      cout << "error: issue with file format" << endl;
      polys.close();
      exit(0);
    }

    // this program is triangles only!!
    if(numVert != 3)
    {
      cout << " This program only makes triangles. Error - too many verts!" << endl;
      polys.close();
      exit(0);
    }

    // loop through to create the triangles
    for(int i = 0; i < numTri; i++)
    {
      //get the next line
      std::getline(polys, line);
      std::istringstream iss2(line);

      // vertices to store
      int index;
      int v1;
      int v2;
      int v3;

      // if it doesnt match this format, there's an error
      if(!(iss2 >> index >> v1 >> v2 >> v3))
      {
        cout << "error: issue with file format" << endl;
        polys.close();
        exit(0);
      }

      // make and save triangle
      vector<VEC2*> triangle;
      triangle.push_back(&_vertices[v1 - 1]);
      triangle.push_back(&_vertices[v2 - 1]);
      triangle.push_back(&_vertices[v3 - 1]);
      _triangles.push_back(TRIANGLE(_material, triangle));
    }
  }
  else // error! shut down
  {
    cout << " Could not open file " << nodeFile << "!!!" << endl;
    nodes.close();
    exit(0);
  }

  polys.close();

  // allocate the state vectors
  _DOFs = 2 * (_vertices.size() - _constrainedVertices.size());

  // if this is not a motion sim with translation
  if (sceneNum)
  {
    size = _unconstrainedVertices.size()*2;
    if(!create_basis)
    {
      basisNoTranslation(filename, cols);
      size = _vertices.size()*2;
    }
  }
  else
  {
    // this function doesn't actually work - currently motion won't be reduced
    setBasisReduction();
    size = _vertices.size()*2;
  }

  // create the mass matrix
  setMassMatrix(!create_basis);

  // set all vectors to zero and determine their size
  VECTOR zeros(size);
  VECTOR z2(_U.cols());

  z2.setZero();
  zeros.setZero();

  _u            = zeros;
  _q            = z2;
  // _ra           = z2;
  // _rv           = z2;
  // _f            = z2;
  _f            = z2;
  _ra           = z2;
  _rv           = z2;
  _fExternal    = zeros;
  _acceleration = zeros;
  _velocity     = zeros;

  // compute the reverse lookup
  computeVertexToIndexTable();

  // if we aren't creating a basis, then create the tensor coefficients for precomputed method
  if(!create_basis)
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
// set new vertex coordinates by adding displacements
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
// set displacements
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

///////////////////////////////////////////////////////////////////////
// set the mass matrix. If reduction, then multiply by U on both sides
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setMassMatrix(bool reduction)
{
  // matrix of size 2Nx2N
  MATRIX M(_vertices.size()*2,_vertices.size()*2);

  //for now we will make every vertex has mass 1
  M.setIdentity();
  M = M*15;

  if (reduction)
  {
    M = M*_U;
    M = _U.transpose() * M;
  }
  _mass = M;
}

///////////////////////////////////////////////////////////////////////
// Create a basis matrix with no translation.
// this is done by reading in a basis file, then performing SVD.
// the basis file is created by doing multiple quasistatic tests
// on the unreduced model, then storing displacements.
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::basisNoTranslation(const char* filename, int basis_cols)
{
  FILE* file = NULL;
  int rows;
  int cols;
  string readin = filename + string(".basis");
  file = fopen(readin.c_str(), "r");
  fscanf(file, "%i %i", &rows, &cols);
  MATRIX U(rows, cols);
  MATRIX svddiag(rows, basis_cols);
  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < cols; j++)
      fscanf(file, "%lf", &U(i,j));
  }

  fclose(file);

  JacobiSVD<MatrixXd> svd( U, ComputeFullV | ComputeFullU );
   svddiag = svd.matrixU();
   svddiag.conservativeResize(svddiag.rows(), basis_cols);
   _U = svddiag;

  // printMatrix(_U);
}

///////////////////////////////////////////////////////////////////////
// Create a basis matrix with translation.
// NOTE: THIS CURRENTLY DOES NOTHING. IGNORE FOR NOW, REWRITE LATER.
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setBasisReduction()
{
  MATRIX U(46,9);
  _U = U;
}

///////////////////////////////////////////////////////////////////////
// Convert reduced coordinates to actual displacement
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::qTou()
{
  _u = _U * _q;
}

///////////////////////////////////////////////////////////////////////
// Add external force to be applied to all vertices
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
// Add force to a single vertex
///////////////////////////////////////////////////////////////////////
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
    if (_restVertices[right][1] > -.85)
      _vertices[right][0] += shear;
  }
}

///////////////////////////////////////////////////////////////////////
// advance the constrained nodes for the stretch test
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::stretch2(const Real stretch)
{
  for (unsigned int x = 0; x < _constrainedVertices.size(); x++)
  {
    int right = _constrainedVertices[x];
    if (_restVertices[right][1] > -.35)
      _vertices[right][1] += stretch;
  }
}

///////////////////////////////////////////////////////////////////////
// Precompute the coefficients for calculating the internal force
// and stiffness matrix.
///////////////////////////////////////////////////////////////////////
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

  // reduction matrix
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
  constTemp.resize(0,0);

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

  // free memory if reduced
  transpose.resize(0,0);
}

///////////////////////////////////////////////////////////////////////
// Compute internal forces using the precomputed method.
// This will find the global force vector of forces on each unrestrained
// vertex.
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeMaterialForces()
{
  //the global vector will be v= [f_0, f_1 ...] where each f_i = [x,y] (column vectors)
  //and forces are only for unconstrained vertices. This means the size of this vector is

  VECTOR displacement = _q;
  // VECTOR displacement = _u;

  VECTOR global_vector(displacement.size());
  global_vector.setZero();

  //add linear term to global vector
  global_vector += (_linearCoef * displacement);

  // compute the cubic term
  TENSOR3 temp = _cubicCoef.modeFourProduct(displacement);
  MATRIX tempMatrix = temp.modeThreeProduct(displacement);

  //add the cubic term to the global vector
  global_vector += (tempMatrix*displacement);

  //compute the quadratic term, add to global vector
  tempMatrix = _cubicquad.modeThreeProduct(displacement);
  global_vector += (tempMatrix * displacement);

  // clear the memory of temporary things
  tempMatrix.resize(0,0);
  temp.clear();

  //store global vector into _f
  _f = global_vector;

}

///////////////////////////////////////////////////////////////////////
// Compute internal forces manually, by looping through triangles.
// This will find the global force vector of forces on each unrestrained
// vertex.
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeUnprecomputedMaterialForces()
{
  //the global vector will be v= [f_0, f_1 ...] where each f_i = [x,y] (column vectors)
  //and forces are only for unconstrained vertices. This means the size of this vector is
  // size(unconstrained vertices)*2 since each has an x,y component
  VECTOR global_vector(_unconstrainedVertices.size()*2);
  global_vector.setZero();

  //loop through triangles, calculate force on each local point
  //convert local vertices to global vertex
  //add force computed into global vertex spot of global_vector
  int n_of_triangles = _triangles.size();

  for(int x = 0; x < n_of_triangles; x++)
  {
    //get our current triangle
    TRIANGLE current = getTriangle(x);

    //find the force vector that's being applied to this triangle
    VECTOR current_force = current.computeForceVector();

    //add the forces into the right place in the global force vector
    for(int y = 0; y < 3; y++)
    {
      //get current triangle vertex y
      VEC2* current_v = _triangles[x].vertex(y);

      //if the vertex we're on is unconstrained, add it to the global vector
      if(_vertexToIndex.find(current_v) != _vertexToIndex.end())
      {
        //find the global index using the vertexToIndex map
        int global_index = _vertexToIndex[current_v];
        global_vector[global_index] += current_force[2*y];
        global_vector[global_index + 1] += current_force[2*y + 1];
      }
    }
  }

  //store global vector into _f
  _f = global_vector;
}

///////////////////////////////////////////////////////////////////////
// A weird collision detection function
///////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////
// A correct collision detection function
///////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////
// Motion step using Euler-Lagrange equation of motion
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::stepMotion(float dt, const VEC2& outerForce, int sceneNum)
{
  //make stiffness Matrix K. size is 2*unrestrained vertices x  2*unrestrained vertices
  // MATRIX K(_u.size(),_u.size() );
  // MATRIX D(_u.size(),_u.size() );

  //reduced
  MATRIX K(_q.size(),_q.size() );
  MATRIX D(_q.size(),_q.size() );

  MATRIX inverse;
  float alpha = 0.01; // constant for damping
  float beta = 0.02;  // constant for damping

  // in barbic, idk if you need to check collision with the wall.
  if (sceneNum == 0)
  {
    printf("checking collision\n");
    checkCollision();
  }

  // Newton Raphson Iteration, but j-max is 1 so no need to write the loop
  //step 1: compute K
  K.setZero();
  computeStiffnessMatrix(K);

  // step 2: compute D
  D = alpha*_mass - beta*K;

  // step 3: calculate f_external
  // VECTOR reducedF = _fExternal;
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
  // _u += dq;
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
  // _velocity = _rv;

  _fExternal.setZero();
  _f.setZero();

}
///////////////////////////////////////////////////////////////////////
// a quasistatic step
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::stepQuasistatic()
{
  //make stiffness Matrix K. size is 2*unrestrained vertices x  2*unrestrained vertices
  MATRIX K(_unconstrainedVertices.size()*2,_unconstrainedVertices.size()*2);
  // MATRIX K(_u.size(),_u.size());
  // MATRIX K(_q.size(),_q.size());

  //step 1: compute K
  K.setZero();
  computeUnprecomputedStiffnessMatrix(K);

  //step 2: compute internal material forces, R(uq)
  computeUnprecomputedMaterialForces();

  // //step 3: external forces transform
  // VECTOR reducedF = _U.transpose() * _fExternal;

  // step 4: form the residual (r = F + E)
  VECTOR r2 = -1*(_f + _fExternal);
  // VECTOR r2 = -1*(_f + reducedF);

  //step 5: compute x = K .inverse().eval()  * r
  // MATRIX inverse2 = K.inverse().eval();
  // VECTOR x2 = inverse2*r2;
  VECTOR x2 = K.colPivHouseholderQr().solve(r2);

  //step 6: add solution x to displamcement vector _u
  _u += x2;
  // _q += x2;
  // qTou();

  //step 7: update all node positions w new displacement vector
  uScatter();

  //reset forces to 0
  _fExternal.setZero();
  _f.setZero();
  K.resize(0,0);

  return true;
}

///////////////////////////////////////////////////////////////////////
// compute the stiffness matrix using precomputed coefficients
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeStiffnessMatrix(MATRIX& K)
{
  //can assume K is the correct size, 2V x 2V

  VECTOR displacement = _q;
  // VECTOR displacement = _u;

  // compute quadratic term and add to K
  TENSOR3 temp = _quadraticCoef.modeFourProduct(displacement);
  K = temp.modeThreeProduct(displacement);

  // compute linear term and add to K
  K += _quadlinear.modeThreeProduct(displacement);

  // add constant term to K
  K += _constCoef;

  // free space
  temp.clear();
}

///////////////////////////////////////////////////////////////////////
// compute the stiffness matrix manually
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeUnprecomputedStiffnessMatrix(MATRIX& K)
{
  //can assume K is the correct size, 2V x 2V

  int n_of_triangles = _triangles.size();

  //go through each triangle and map vertices to their global vertices. if
  //unrestrained, then add to the poisition in the global K.
  for(int x = 0; x < n_of_triangles; x++)
  {
    TRIANGLE current = getTriangle(x);

    //find the force vector that's being applied to this triangle
    MATRIX current_force = current.computeForceJacobian();

    for(int i = 0; i < 3; i++)
    {
      //get current triangle vertex y
      VEC2* current_iv = _triangles[x].vertex(i);

      //if the vertex we're on is unconstrained, add it to the global K
      if(_vertexToIndex.find(current_iv) != _vertexToIndex.end())
      {
        //find the global index using the vertexToIndex map
        int global_index_i = _vertexToIndex[current_iv];

        for(int j = 0; j < 3; j++)
        {
          VEC2* current_jv = _triangles[x].vertex(j);

          if(_vertexToIndex.find(current_jv) != _vertexToIndex.end())
          {
            int global_index_j = _vertexToIndex[current_jv];
            K(global_index_i, global_index_j) += current_force(2*i, 2*j);
            K(global_index_i, global_index_j + 1) += current_force(2*i, 2*j + 1);
            K(global_index_i + 1, global_index_j) += current_force(2*i + 1, 2*j);
            K(global_index_i + 1, global_index_j + 1) += current_force(2*i + 1, 2*j + 1);
          }
        }
      }
    }
  }
}
