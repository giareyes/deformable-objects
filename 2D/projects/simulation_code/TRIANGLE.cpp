#include "TRIANGLE.h"

#include <iostream>
using namespace std;

TRIANGLE::TRIANGLE(MATERIAL* material, const vector<VEC2*>& vertices) :
  _material(material)
{
  assert(vertices.size() == 3);
  _vertices[0] = vertices[0];
  _vertices[1] = vertices[1];
  _vertices[2] = vertices[2];

  // store these as the rest pose, for now
  for (unsigned int x = 0; x < 3; x++)
    _restPose[x] = *_vertices[x];

    VEC2 e0 = _restPose[1] - _restPose[0];
    VEC2 e1 = _restPose[2] - _restPose[0];

    // put e0, e1 into Dm
    _Dm.col(0) = e0;
    _Dm.col(1) = e1;

    _lambda = _material->getLambda();
    _mu = _material->getMu();

    // create the constant coefficient matrix below:
    // this is equal to (pfpu)T *(mu*I - (lambda + mu)d(detF)/dF)* (pfpu)
    _pfpu = pFpuVectorized();
    MATRIX constMatrix(4,4);
    constMatrix.setIdentity();
    constMatrix = constMatrix*_mu;
    constMatrix(0,3) = -1*(_lambda + _mu);
    constMatrix(1,2) = (_lambda + _mu);
    constMatrix(2,1) = (_lambda + _mu);
    constMatrix(3,0) = -1*(_lambda + _mu);

    _constCoef = _pfpu.transpose() * constMatrix * _pfpu;
    constMatrix.resize(0,0);

    // multiply by -1*restArea
    _constCoef = -1*restArea()*_constCoef;

    // create the _quadraticCoef for the stiffness matrix polynomial
    // it must be multiplied in all dimensions by (pfpu)^T
    TENSOR4 quad(4,4,4,4);
    // column 0
    // row 0
    quad._tensor[0]._tensor[0](3,3) = 1;
    // row 1
    quad._tensor[0]._tensor[1](3,2) = -1;
    // row 2
    quad._tensor[0]._tensor[2](3,1) = -1;
    // row 3
    quad._tensor[0]._tensor[3](2,1) = -1;
    quad._tensor[0]._tensor[3](3,0) = 2;

    // column 1
    // row 0
    quad._tensor[1]._tensor[0](2,3) = -1;
    // row 1
    quad._tensor[1]._tensor[1](2,2) = 1;
    // row 2
    quad._tensor[1]._tensor[2](2,1) = 2;
    quad._tensor[1]._tensor[2](3,0) = -1;
    // row 3
    quad._tensor[1]._tensor[3](2,0) = -1;

    // column 2
    // row 0
    quad._tensor[2]._tensor[0](1,3) = -1;
    // row 1
    quad._tensor[2]._tensor[1](0,3) = -1;
    quad._tensor[2]._tensor[1](1,2) = 2;
    // row 2
    quad._tensor[2]._tensor[2](1,1) = 1;
    // row 3
    quad._tensor[2]._tensor[3](1,0) = -1;

    // column 3
    // row 0
    quad._tensor[3]._tensor[0](0,3) = 2;
    quad._tensor[3]._tensor[0](1,2) = -1;
    // row 1
    quad._tensor[3]._tensor[1](0,2) = -1;
    // row 2
    quad._tensor[3]._tensor[2](0,1) = -1;
    // row 3
    quad._tensor[3]._tensor[3](0,0) = 1;

    _quadraticCoef = quad;
    _quadraticCoef *= (-1*restArea()*_lambda);

    // multiply in all dimensions by (pfpu)^T
    MATRIX pfputrans = _pfpu.transpose();
    _quadraticCoef = _quadraticCoef.modeFourProduct(pfputrans);
    _quadraticCoef = _quadraticCoef.modeThreeProduct(pfputrans);
    _quadraticCoef = _quadraticCoef.modeTwoProduct(pfputrans);
    _quadraticCoef = _quadraticCoef.modeOneProduct(pfputrans);

    // create the cubic term coefficient for the internal force polynomial.
    // again, all dimensions must be multiplied by (pfpu)^T
    TENSOR4 cubic(4,4,4,4);

    for(int i = 0; i < 4; i++)
    {
      cubic._tensor[i]._tensor[i](3-i, 3-i) = 1;
    }

    cubic._tensor[0]._tensor[1](3, 2) = -1;
    cubic._tensor[1]._tensor[0](2, 3) = -1;

    cubic._tensor[2]._tensor[3](1, 0) = -1;
    cubic._tensor[3]._tensor[2](0, 1) = -1;

    _cubicCoef = cubic;
    _cubicCoef *= (-1*restArea()*_lambda);

    // multiply in all dimensions by (pfpu)^T
    _cubicCoef = _cubicCoef.modeFourProduct(pfputrans);
    _cubicCoef = _cubicCoef.modeThreeProduct(pfputrans);
    _cubicCoef = _cubicCoef.modeTwoProduct(pfputrans);
    _cubicCoef = _cubicCoef.modeOneProduct(pfputrans);

    // we no longer need this, so free memeory!
    pfputrans.resize(0,0);

    // add on the terms that are created when we shift our equations to handle displacement instead of position

    // first we must get all rest positions and store it in a vector of size 6
    VECTOR pos(6);
    for(int i = 0; i < 3; i++)
    {
      pos[2*i] = (_restPose[i])[0];
      pos[2*i + 1] = (_restPose[i])[1];
    }

    //create an initial quadratic term for the internal force polynomial
    _cubic2 = _cubicCoef.modeFourProduct(pos); //D x_4 x;

    // this initial quadratic term is the base for the linear coefficient of the internal force polynomial as well
    _linearCoef = _cubic2.modeThreeProduct(pos); // D x_4 x x_3 x
    _linearCoef += _cubic2.modeTwoProduct(pos); // D x_4 x x_2 x

    // create the next quadratic term for the internal force polynomial
    TENSOR3 tempCubic = _cubicCoef.modeThreeProduct(pos); // D x_3 x

    // this is also used as a base for one of the linear terms
    _linearCoef += tempCubic.modeTwoProduct(pos); // D x_3 x x_2 x

    // now add the constant coefficient from the stiffness matrix to the linear coefficient
    _linearCoef += _constCoef;

    // add the last quadratic term to our constant
    _cubic2 += tempCubic;
    _cubic2 += _cubicCoef.modeTwoProduct(pos);

    // we are done with tempCubic, so free the memory!
    tempCubic.clear();

    // create our linear term for the stiffness matrix polynomial
    _quad2 = _quadraticCoef.modeFourProduct(pos); // C x_4 x;

    // this linear term is used as a base for what's added to the constant coefficient
    _constCoef += _quad2.modeThreeProduct(pos); // C x_4 x x_3 x;

    // add the last linear term to the stiffness matrix coefficent
    _quad2 += _quadraticCoef.modeThreeProduct(pos); // C x_3 x;


}

///////////////////////////////////////////////////////////////////////
// get the deformation gradient at Xi
///////////////////////////////////////////////////////////////////////
MATRIX2 TRIANGLE::computeF() const
{
  MATRIX2 F; //force gradient; will be returned
  MATRIX2 Ds; // spatial matrix
  MATRIX2 DmInverse;

  //deformed e0 and e1
  VEC2 de0 = *_vertices[1] - *_vertices[0];
  VEC2 de1 = *_vertices[2] - *_vertices[0];

  //put de0 and de1 into our spacial matrix
  Ds.col(0) = de0;
  Ds.col(1) = de1;

  //multiply Ds and Dm inverse to calcuate F
  DmInverse = _Dm.inverse().eval();
  F = Ds*DmInverse;

  return F;
}

///////////////////////////////////////////////////////////////////////
// take the average of the vertices
///////////////////////////////////////////////////////////////////////
VEC2 TRIANGLE::vertexAverage()
{
  VEC2 final = (*_vertices[0]);
  for (int x = 1; x < 3; x++)
    final += (*_vertices[x]);

  return final * 1.0 / 3.0;
}

///////////////////////////////////////////////////////////////////////
// use a precomputed Linear Coefficient
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::precomputedLinearCoef()
{
  // figure out how to get lamda and mu
  VECTOR linearForce(6);
  VECTOR pos(6);
  for(int i = 0; i < 3; i++)
  {
    pos[2*i] = (*_vertices[i])[0];
    pos[2*i + 1] = (*_vertices[i])[1];
  }

  linearForce = _constCoef*pos;
  return linearForce;
}

VECTOR TRIANGLE::precomputedCubicCoef() {
  // figure out how to get lamda and mu
  VECTOR cubicForce(6);
  VECTOR pos(6);
  for(int i = 0; i < 3; i++)
  {
    pos[2*i] = (*_vertices[i])[0];
    pos[2*i + 1] = (*_vertices[i])[1];
  }

  TENSOR3 cubicTensor = _cubicCoef.modeFourProduct(pos);
  MATRIX cubicMatrix = cubicTensor.modeThreeProduct(pos);
  cubicForce = cubicMatrix * pos;

  cubicTensor.clear();
  cubicMatrix.resize(0,0);

  return cubicForce;
}

MATRIX TRIANGLE::precomputedQuadCoef()
{
  MATRIX quadForce(6,6);

  VECTOR pos(6);
  for(int i = 0; i < 3; i++)
  {
    pos[2*i] = (*_vertices[i])[0];
    pos[2*i + 1] = (*_vertices[i])[1];
  }

  TENSOR3 quadTensor = _quadraticCoef.modeFourProduct(pos);

  quadForce = quadTensor.modeThreeProduct(pos);

  quadTensor.clear();

  return  quadForce;
}

///////////////////////////////////////////////////////////////////////
// populate the force vector
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::computeForceVector()
{
  VECTOR forceVector(6);
  MATRIX2 F = computeF();
  MATRIX2 pk1 = _material->PK1(F);
  MATRIX dfdu(4,6);

//multiply pk1 by -1 and area as discussed in OH
  pk1 = -1*restArea()*pk1;

  //get vectorized df/du
  dfdu = pFpuVectorized();

  forceVector = dfdu.transpose()*vectorize(pk1);

  return forceVector;
}

///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::computeForceJacobian()
{
  MATRIX2 F = computeF();
  MATRIX jacobian(6,6);
  MATRIX dfdu(4,6);

  //get 4x4 matrix of derivative of pk1. This is matrix B from class
  MATRIX dpdf = _material->DPDF(F);

  //multiply dpdf by -1 and area as discussed in OH
  dpdf = -1*restArea()*dpdf;

  //get vectorized df/du. This is matrix A from class
  dfdu = pFpuVectorized();

  // d^2 psi/ dx^2 = A transpose * B * A
  jacobian = dfdu.transpose() * dpdf * dfdu;

  return jacobian;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
vector<MATRIX> TRIANGLE::pFpu()
{
  vector<MATRIX> result; //the end 3rd order tensor
  MATRIX2 DmInverse;
  DmInverse = _Dm.inverse().eval();

  //compute dF/x0 of [[[ x2 - x0], [x3 - x1]], [[ x4 - x0], [x5 - x1]]]
  MATRIX2 x0;
  x0.setZero();
  x0(0,0) = -1;
  x0(0,1) = -1;
  x0 = x0*DmInverse;
  result.push_back(x0);

  //compute dF/x1 of [[[ x2 - x0], [x3 - x1]], [[ x4 - x0], [x5 - x1]]]
  MATRIX2 x1;
  x1.setZero();
  x1(1,0) = -1;
  x1(1,1) = -1;
  x1 = x1*DmInverse;
  result.push_back(x1);

  //compute dF/x2 of [[[ x2 - x0], [x3 - x1]], [[ x4 - x0], [x5 - x1]]]
  MATRIX2 x2;
  x2.setZero();
  x2(0,0) = 1;
  x2 = x2*DmInverse;
  result.push_back(x2);

  //compute dF/x3 of [[[ x2 - x0], [x3 - x1]], [[ x4 - x0], [x5 - x1]]]
  MATRIX2 x3;
  x3.setZero();
  x3(1,0) = 1;
  x3 = x3*DmInverse;
  result.push_back(x3);

  //compute dF/x4 of [[[ x2 - x0], [x3 - x1]], [[ x4 - x0], [x5 - x1]]]
  MATRIX2 x4;
  x4.setZero();
  x4(0,1) = 1;
  x4 = x4*DmInverse;
  result.push_back(x4);

  //compute dF/x5 of [[[ x2 - x0], [x3 - x1]], [[ x4 - x0], [x5 - x1]]]
  MATRIX2 x5;
  x5.setZero();
  x5(1,1) = 1;
  x5 = x5*DmInverse;
  result.push_back(x5);

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::pFpuVectorized()
{
  //take result from pFpu and vectorize it using function from EXTRAFUNCTIONS.cpp
  vector<MATRIX> pfpu = pFpu();
  MATRIX vectorized(4,6);
  for(int x = 0; x < 6; x++)
  {
    vectorized.col(x) = vectorize(pfpu[x]);
  }
  return vectorized;
}

///////////////////////////////////////////////////////////////////////
// compute rest area of this triangle
///////////////////////////////////////////////////////////////////////
Real TRIANGLE::restArea() const
{
  // get triangle normal in R^3
  VEC3 triangleNormal;
  VEC3 restPose3[3];
  for (int x = 0; x < 3; x++)
  {
    restPose3[x].setZero();
    for (int y = 0; y < 2; y++)
      restPose3[x][y] = _restPose[x][y];
  }
  triangleNormal = (restPose3[2] - restPose3[0]).cross(restPose3[1] - restPose3[0]);
  return triangleNormal.norm() * 0.5;
}

///////////////////////////////////////////////////////////////////////
// compute deformed area of this triangle
///////////////////////////////////////////////////////////////////////
Real TRIANGLE::area() const
{
  // get triangle normal in R^3
  VEC3 triangleNormal;
  VEC3 pose3[3];
  for (int x = 0; x < 3; x++)
  {
    pose3[x].setZero();
    for (int y = 0; y < 2; y++)
      pose3[x][y] = (*_vertices[x])[y];
  }
  triangleNormal = (pose3[2] - pose3[0]).cross(pose3[1] - pose3[0]);
  return triangleNormal.norm() * 0.5;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::getDisplacements() const
{
  VECTOR result(6);
  result[0] = (*_vertices[0])[0] - _restPose[0][0];
  result[1] = (*_vertices[0])[1] - _restPose[0][1];
  result[2] = (*_vertices[1])[0] - _restPose[1][0];
  result[3] = (*_vertices[1])[1] - _restPose[1][1];
  result[4] = (*_vertices[2])[0] - _restPose[2][0];
  result[5] = (*_vertices[2])[1] - _restPose[2][1];

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE::setDisplacements(const VECTOR& u)
{
  (*_vertices[0])[0] = _restPose[0][0] + u[0];
  (*_vertices[0])[1] = _restPose[0][1] + u[1];
  (*_vertices[1])[0] = _restPose[1][0] + u[2];
  (*_vertices[1])[1] = _restPose[1][1] + u[3];
  (*_vertices[2])[0] = _restPose[2][0] + u[4];
  (*_vertices[2])[1] = _restPose[2][1] + u[5];
}
