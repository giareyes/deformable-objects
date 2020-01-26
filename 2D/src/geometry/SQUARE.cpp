#include "SQUARE.h"
#include <GLUT/glut.h>
#include "TIMER.h"
#include "HYPER_TAN.h"
#include "INVERTIBLE.h"

#include <iostream>
using namespace std;

SQUARE::SQUARE(MATERIAL* material, const vector<VEC2*>& vertices) :
  _material(material)
{
  assert(vertices.size() == 4);
  _vertices[0] = vertices[0];
  _vertices[1] = vertices[1];
  _vertices[2] = vertices[2];
  _vertices[3] = vertices[3];

  //_vertices[0] = VEC2(-1, -1);
  //_vertices[1] = VEC2(1, -1);
  //_vertices[2] = VEC2(1, 1);
  //_vertices[3] = VEC2(-1, 1);
  _Dm = vertexMatrix();
  
  // pre-cache some inverses for F
  const Real weight = 1.0;
  const Real invSqrt3 = 1.0 / sqrt(3.0);
  const VEC2 Xi[] = {invSqrt3 * VEC2(-1,-1),
                     invSqrt3 * VEC2(1,-1),
                     invSqrt3 * VEC2(1,1),
                     invSqrt3 * VEC2(-1,1)};
  for (int x = 0; x < 4; x++)
  {
    _H[x]       = computeH(Xi[x]);
    _DmH_Inv[x] = (_Dm * _H[x]).inverse();

    MATRIX DmHg = _Dm * _H[x];
    Real dXdXi = DmHg.determinant();
    MATRIX DmHgInverse = DmHg.inverse();
    _Bmg[x] = DmHgInverse.transpose() * _H[x].transpose() * dXdXi * weight;

    _DFDuA[x] = _H[x] * (_Dm * _H[x]).inverse();
  }

  // precache DFDu
  for (int g = 0; g < 4; g++)
    for (int i = 0; i < 8; i++)
    {
      MATRIX DFDu_i = computeDFDu(g, i);
      _flatDFDu[g][i] = flatten(DFDu_i);
    }

  for (int x = 0; x < 4; x++)
    _restPose[x] = *_vertices[x];

  _forces.resize(2,4);

  // init the inversion variables
  for (int x = 0; x < 4; x++)
  {
    _oldClamped[x] = false;
    _oldDirection[x].setZero();
    _newClamped[x] = false;
    _newDirection[x].setZero();
  }
}

///////////////////////////////////////////////////////////////////////
// Draw the square to GL
///////////////////////////////////////////////////////////////////////
void SQUARE::draw()
{
  /*
  glPointSize(10.0);
  glColor4f(0,0,1,1);
  glBegin(GL_POINTS);
    for (int x = 0; x < 4; x++)
      glVertex2f((*_vertices[x])[0], (*_vertices[x])[1]);
  glEnd();
  */
 
  // draw a filled in quad
  //glColor4f(0.5, 0.5, 0.5, 1);
  //glColor4f(0.5, 1.0, 0.5, 0.5);
  
  //VEC3 colormap = ramp(1 - _conditionNumber);
  VEC3 colormap = ramp(_conditionNumber);
  glColor4f(colormap[0], colormap[1], colormap[2], 1.0);
  glBegin(GL_TRIANGLES);
    glVertex2f((*_vertices[0])[0], (*_vertices[0])[1]);
    glVertex2f((*_vertices[1])[0], (*_vertices[1])[1]);
    glVertex2f((*_vertices[2])[0], (*_vertices[2])[1]);

    glVertex2f((*_vertices[2])[0], (*_vertices[2])[1]);
    glVertex2f((*_vertices[3])[0], (*_vertices[3])[1]);
    glVertex2f((*_vertices[0])[0], (*_vertices[0])[1]);
  glEnd(); 

  glLineWidth(1);

  // draw the outline
  glColor4f(0,0,0,1);
  glBegin(GL_LINE_STRIP);
    for (int x = 0; x < 4; x++)
      glVertex2f((*_vertices[x])[0], (*_vertices[x])[1]);
    glVertex2f((*_vertices[0])[0], (*_vertices[0])[1]);
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// Draw the inversion direction to GL
///////////////////////////////////////////////////////////////////////
void SQUARE::drawInversionDirection()
{
  VEC2 center = *_vertices[0];
  for (int x = 1; x < 4; x++)
    center += *_vertices[x];
  center *= 0.25;

  VEC2 end = center + _drawDirection[0];

  glLineWidth(2);
  glBegin(GL_LINES);
    glColor4f(0,1,0,1);
    glVertex2f(center[0], center[1]);
    glVertex2f(end[0], end[1]);
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// Draw the square forces to GL
///////////////////////////////////////////////////////////////////////
void SQUARE::drawForces()
{
  if (_forceCache.size() != 8)
    return;

  glLineWidth(2);
  glBegin(GL_LINES);
    for (int x = 0; x < 4; x++)
    {
      glColor4f(1,0,0,1);
      glVertex2f((*_vertices[x])[0], (*_vertices[x])[1]);

      Real dt = 0.001;
      VEC2 force(_forceCache[2 * x], _forceCache[2 * x + 1]);

      force.normalize();
      //VEC2 forceEnd = (*_vertices[x]) + _forces.col(x) * dt;
      VEC2 forceEnd = *_vertices[x] + force;
      
      glColor4f(1,0,0,0);
      glVertex2f(forceEnd[0], forceEnd[1]);
    }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// get the shape function derivative matrix H for computing F
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeH(const VEC2& Xi) const
{
  MATRIX final(4,2);

  final(0,0) = -(1.0 - Xi[1]);
  final(1,0) =  (1.0 - Xi[1]);
  final(2,0) =  (1.0 + Xi[1]);
  final(3,0) = -(1.0 + Xi[1]);
  
  final(0,1) = -(1.0 - Xi[0]);
  final(1,1) = -(1.0 + Xi[0]);
  final(2,1) =  (1.0 + Xi[0]);
  final(3,1) =  (1.0 - Xi[0]);

  return final * 0.25;
}

///////////////////////////////////////////////////////////////////////
// build the material coordinate matrix
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::vertexMatrix() const
{
  MATRIX final(2,4);
  final(0,0) = (*_vertices[0])[0];
  final(0,1) = (*_vertices[1])[0];
  final(0,2) = (*_vertices[2])[0];
  final(0,3) = (*_vertices[3])[0];
  
  final(1,0) = (*_vertices[0])[1];
  final(1,1) = (*_vertices[1])[1];
  final(1,2) = (*_vertices[2])[1];
  final(1,3) = (*_vertices[3])[1];

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the deformation gradient at Xi, the normalized position
// inside a [-1,-1], [1,1] square
///////////////////////////////////////////////////////////////////////
MATRIX2 SQUARE::computeF(const VEC2& Xi) const
{
  MATRIX H = computeH(Xi);
  MATRIX Ds = vertexMatrix();

  MATRIX2 DmH_Inv = (_Dm * H).inverse();

  return (Ds * H) * DmH_Inv;
}

///////////////////////////////////////////////////////////////////////
// scale whole square
///////////////////////////////////////////////////////////////////////
void SQUARE::scale(const VEC2& scalar)
{
  for (int x = 0; x < 4; x++)
  {
    (*_vertices[x])[0] *= scalar[0];
    (*_vertices[x])[1] *= scalar[1];
  }
}

///////////////////////////////////////////////////////////////////////
// translate the whole square
///////////////////////////////////////////////////////////////////////
void SQUARE::translate(const VEC2& translation)
{
  for (int x = 0; x < 4; x++)
    (*_vertices[x]) += translation;
}

///////////////////////////////////////////////////////////////////////
// rotate the whole square
///////////////////////////////////////////////////////////////////////
void SQUARE::rotate(const Real theta)
{
  VEC2 center = vertexAverage();

  MATRIX2 rotation;
  rotation(0,0) = rotation(1,1) = cos(theta);
  rotation(0,1) = -sin(theta);
  rotation(1,0) = sin(theta);

  for (int x = 0; x < 4; x++)
    *_vertices[x] = rotation * ((*_vertices[x]) - center) + center;
}

///////////////////////////////////////////////////////////////////////
// take the average of the vertices
///////////////////////////////////////////////////////////////////////
VEC2 SQUARE::vertexAverage()
{
  VEC2 final = *_vertices[0];
  for (int x = 1; x < 4; x++)
    final += *_vertices[x];

  return final * 0.25;
}

///////////////////////////////////////////////////////////////////////
// get the average F at all the Gauss points
///////////////////////////////////////////////////////////////////////
MATRIX2 SQUARE::gaussF()
{
  const Real invSqrt3 = 1.0 / sqrt(3.0);
  /*
  const VEC2 gaussPoints[] = {invSqrt3 * (*_vertices[0]),
                              invSqrt3 * (*_vertices[1]),
                              invSqrt3 * (*_vertices[2]),
                              invSqrt3 * (*_vertices[3])};
                              */
  const VEC2 gaussPoints[] = {invSqrt3 * VEC2(-1,-1),
                              invSqrt3 * VEC2(1,-1),
                              invSqrt3 * VEC2(1,1),
                              invSqrt3 * VEC2(-1,1)};

  MATRIX2 final = computeF(gaussPoints[0]);
  for (int x = 1; x < 4; x++)
    final += computeF(gaussPoints[x]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// populate the force vector
///////////////////////////////////////////////////////////////////////
void SQUARE::computeForces()
{
  // build the Gauss points in the canonical volume
  const Real invSqrt3 = 1.0 / sqrt(3.0);
  const VEC2 gaussPoints[] = {invSqrt3 * VEC2(-1,-1),
                              invSqrt3 * VEC2(1,-1),
                              invSqrt3 * VEC2(1,1),
                              invSqrt3 * VEC2(-1,1)};

  _forces.setZero();
  for (int x = 0; x < 4; x++)
    _forces += computeForces(gaussPoints[x]);
}

///////////////////////////////////////////////////////////////////////
// populate the force vector
///////////////////////////////////////////////////////////////////////
VECTOR SQUARE::computeForceVector()
{
  // build the Gauss points in the canonical volume
  const Real invSqrt3 = 1.0 / sqrt(3.0);
  const VEC2 gaussPoints[] = {invSqrt3 * VEC2(-1,-1),
                              invSqrt3 * VEC2(1,-1),
                              invSqrt3 * VEC2(1,1),
                              invSqrt3 * VEC2(-1,1)};

  VECTOR forceVector(8);
  forceVector.setZero();
  for (int x = 0; x < 4; x++)
    forceVector += computeForceVector(gaussPoints[x]);

  return forceVector;
}

///////////////////////////////////////////////////////////////////////
// populate the force vector
///////////////////////////////////////////////////////////////////////
VECTOR SQUARE::computeForceVectorFast()
{
  VECTOR forceVector(8);
  forceVector.setZero();
  for (int x = 0; x < 4; x++)
    forceVector += computeForceVectorFast(x);

  _forceCache = forceVector;

  return forceVector;
}

///////////////////////////////////////////////////////////////////////
// populate the force vector
///////////////////////////////////////////////////////////////////////
VECTOR SQUARE::computeInvertibleForceVector()
{
  VECTOR forceVector(8);
  forceVector.setZero();
  for (int x = 0; x < 4; x++)
    forceVector += computeInvertibleForceVector(x);

  _forceCache = forceVector;

  return forceVector;
}

///////////////////////////////////////////////////////////////////////
// get the forces at a Gauss point
///////////////////////////////////////////////////////////////////////
VECTOR SQUARE::computeForceVector(const VEC2& gaussPoint)
{
  // the Gauss weight is one for the 2x2 scheme
  Real weight = 1.0;

  // all of Bmg can be precomputed, but that is for later
  MATRIX H = computeH(gaussPoint);
  MATRIX DmHg = _Dm * H;
  Real dXdXi = DmHg.determinant();

  MATRIX DmHgInverse = DmHg.inverse();
  MATRIX Bmg = DmHgInverse.transpose() * H.transpose() * dXdXi * weight;

  VECTOR vectorBmg = flatten(Bmg);

  // get the first PK of the actual material
  //MATRIX P = firstPiolaKirchhoffStVK(gaussPoint);
  MATRIX2 F = computeF(gaussPoint);
  MATRIX P = _material->PK1(F);

  MATRIX copiedP(8,8);
  copiedP.setZero();

  for (int x = 0; x < 2; x++)
    for (int y = 0; y < 2; y++)
    {
      copiedP(x,y) = P(x,y);
      copiedP(2 + x, 2 + y) = P(x,y);
      copiedP(4 + x, 4 + y) = P(x,y);
      copiedP(6 + x, 6 + y) = P(x,y);
    }

  return (-1.0) * copiedP * vectorBmg;
}

///////////////////////////////////////////////////////////////////////
// get the forces at a Gauss point
///////////////////////////////////////////////////////////////////////
VECTOR SQUARE::computeForceVectorFast(const int gaussPoint)
{
  /*
  // all of Bmg can be precomputed, but that is for later
  MATRIX H = computeH(gaussPoint);
  MATRIX DmHg = _Dm * H;
  Real dXdXi = DmHg.determinant();
  MATRIX DmHgInverse = DmHg.inverse();
  MATRIX Bmg = DmHgInverse.transpose() * H.transpose() * dXdXi * weight;
  */
  
  MATRIX& H = _H[gaussPoint];
  MATRIX& Bmg = _Bmg[gaussPoint];
  VECTOR vectorBmg = flatten(Bmg);

  // get the first PK of the actual material
  //MATRIX P = firstPiolaKirchhoffStVK(gaussPoint);
  MATRIX Ds = vertexMatrix();
  MATRIX F = (Ds * H) * _DmH_Inv[gaussPoint];

  MATRIX P = _material->PK1(F);

  MATRIX copiedP(8,8);
  copiedP.setZero();

  for (int x = 0; x < 2; x++)
    for (int y = 0; y < 2; y++)
    {
      copiedP(x,y) = P(x,y);
      copiedP(2 + x, 2 + y) = P(x,y);
      copiedP(4 + x, 4 + y) = P(x,y);
      copiedP(6 + x, 6 + y) = P(x,y);
    }

  return (-1.0) * copiedP * vectorBmg;
}

///////////////////////////////////////////////////////////////////////
// get the forces at a Gauss point
///////////////////////////////////////////////////////////////////////
VECTOR SQUARE::computeInvertibleForceVector(const int gaussPoint)
{
  MATRIX& H = _H[gaussPoint];
  MATRIX& Bmg = _Bmg[gaussPoint];
  VECTOR vectorBmg = flatten(Bmg);

  // get the first PK of the actual material
  MATRIX Ds = vertexMatrix();
  MATRIX F = (Ds * H) * _DmH_Inv[gaussPoint];

  // get the invertible material. The pointer IS invertible, right?
  // this will throw an ugly exception if that is not true.
  INVERTIBLE* invertible = dynamic_cast<INVERTIBLE*>(_material);

  // use inversion consistency
  MATRIX P = invertible->PK1(F, _oldDirection[gaussPoint], _oldClamped[gaussPoint],
                                _newDirection[gaussPoint], _newClamped[gaussPoint]);

  MATRIX copiedP(8,8);
  copiedP.setZero();

  for (int x = 0; x < 2; x++)
    for (int y = 0; y < 2; y++)
    {
      copiedP(x,y) = P(x,y);
      copiedP(2 + x, 2 + y) = P(x,y);
      copiedP(4 + x, 4 + y) = P(x,y);
      copiedP(6 + x, 6 + y) = P(x,y);
    }

  return (-1.0) * copiedP * vectorBmg;
}

///////////////////////////////////////////////////////////////////////
// get the forces at a Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeForces(const VEC2& gaussPoint)
{
  // the Gauss weight is one for the 2x2 scheme
  Real weight = 1.0;

  // all of Bmg can be precomputed, but that is for later
  MATRIX H = computeH(gaussPoint);
  MATRIX DmHg = _Dm * H;
  Real dXdXi = DmHg.determinant();
  MATRIX DmHgInverse = DmHg.inverse();
  MATRIX Bmg = DmHgInverse.transpose() * H.transpose() * dXdXi * weight;

  // get the first PK of the actual material
  //MATRIX P = firstPiolaKirchhoffStVK(gaussPoint);
  //MATRIX P = firstPiolaKirchhoffCorotational(gaussPoint);
  MATRIX2 F = computeF(gaussPoint);
  MATRIX P = _material->PK1(F);

  return (-1.0) * P * Bmg;
}

///////////////////////////////////////////////////////////////////////
// get the gradient of PK2 at a Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeDPDu()
{
  // build the Gauss points in the canonical volume
  const Real invSqrt3 = 1.0 / sqrt(3.0);
  const VEC2 gaussPoints[] = {invSqrt3 * VEC2(-1,-1),
                              invSqrt3 * VEC2(1,-1),
                              invSqrt3 * VEC2(1,1),
                              invSqrt3 * VEC2(-1,1)};
  MATRIX final(4,8);
  final.setZero();

  for (int i = 0; i < 8; i++)
  {
    MATRIX col = computeDPDu(gaussPoints[0],i);
    for (int x = 1; x < 4; x++)
      col += computeDPDu(gaussPoints[x],i);
    final.col(i) = SQUARE::flatten(col);
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the forces Jacobian at a Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeForceJacobian()
{
  TIMER functionTimer(__FUNCTION__);
  // build the Gauss points in the canonical volume
  const Real invSqrt3 = 1.0 / sqrt(3.0);
  const VEC2 gaussPoints[] = {invSqrt3 * VEC2(-1,-1),
                              invSqrt3 * VEC2(1,-1),
                              invSqrt3 * VEC2(1,1),
                              invSqrt3 * VEC2(-1,1)};

  MATRIX final;
  final = computeForceJacobian(gaussPoints[0]);
  for (int x = 1; x < 4; x++)
    final += computeForceJacobian(gaussPoints[x]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the forces Jacobian at a Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeForceJacobianFast()
{
  TIMER functionTimer(__FUNCTION__);
  /*
  // build the Gauss points in the canonical volume
  const Real invSqrt3 = 1.0 / sqrt(3.0);
  const VEC2 gaussPoints[] = {invSqrt3 * VEC2(-1,-1),
                              invSqrt3 * VEC2(1,-1),
                              invSqrt3 * VEC2(1,1),
                              invSqrt3 * VEC2(-1,1)};
                              */

  MATRIX final;
  //final = computeForceJacobian(gaussPoints[0]);
  final = computeForceJacobianFast(0);
  for (int x = 1; x < 4; x++)
    final += computeForceJacobianFast(x);

  //EigenSolver<MATRIX> eigensolver(final);
  //cout << " K eigenvalues: " << -1 * eigensolver.eigenvalues().real().transpose() << endl;
  //cout << " imaginary magnitude: " << eigensolver.eigenvalues().imag().norm() << endl;
#if 0
  // cache the condition number
  SelfAdjointEigenSolver<MATRIX> eigensolver(final);
  //JacobiSVD<MATRIX> eigensolver(final);
  MATRIX vectors = eigensolver.eigenvectors();

  /*
  _minEig = fabs(eigensolver.eigenvalues()[0]);
  _maxEig = fabs(eigensolver.eigenvalues()[0]);
  for (int i = 1; i < 8; i++)
  {
    if (fabs(eigensolver.eigenvalues()[i]) < _minEig)
      _minEig = fabs(eigensolver.eigenvalues()[i]);
    if (fabs(eigensolver.eigenvalues()[i]) > _maxEig)
      _maxEig = fabs(eigensolver.eigenvalues()[i]);
  }

  if (_minEig < 1e-16)
    _minEig = 1e-16;
    */

  _minEig = (eigensolver.eigenvalues()[0]);
  _maxEig = (eigensolver.eigenvalues()[0]);
  for (int i = 1; i < 8; i++)
  {
    if ((eigensolver.eigenvalues()[i]) < _minEig)
      _minEig = (eigensolver.eigenvalues()[i]);
    if ((eigensolver.eigenvalues()[i]) > _maxEig)
      _maxEig = (eigensolver.eigenvalues()[i]);
  }

  //if (fabs(_minEig) < 1e-16)
  //  _minEig = 1e-16;
  //_conditionNumber = fabs(_maxEig / _minEig);
  //_conditionNumber = fabs(_maxEig);
  //_conditionNumber = _maxEig;
  _conditionNumber = psi();
#endif

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the forces Jacobian at a Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeInvertibleForceJacobian()
{
  TIMER functionTimer(__FUNCTION__);

  MATRIX final;
  final = computeInvertibleForceJacobian(0);
  for (int x = 1; x < 4; x++)
    final += computeInvertibleForceJacobian(x);

  return final;
}

///////////////////////////////////////////////////////////////////////
// get PK2 over all the Gauss points
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::secondPiolaKirchhoff()
{
  // build the Gauss points in the canonical volume
  const Real invSqrt3 = 1.0 / sqrt(3.0);
  const VEC2 gaussPoints[] = {invSqrt3 * VEC2(-1,-1),
                              invSqrt3 * VEC2(1,-1),
                              invSqrt3 * VEC2(1,1),
                              invSqrt3 * VEC2(-1,1)};

  MATRIX final;
  MATRIX2 F = computeF(gaussPoints[0]);
  final = _material->PK2(F);
  //final = secondPiolaKirchhoffStVK(gaussPoints[0]);
  for (int x = 1; x < 4; x++)
  {
    MATRIX2 F = computeF(gaussPoints[x]);
    final += _material->PK2(F);
    //final += secondPiolaKirchhoffStVK(gaussPoints[x]);
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the forces Jacobian at a Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeForceJacobian(const VEC2& gaussPoint)
{
  // the Gauss weight is one for the 2x2 scheme
  Real weight = 1.0;

  // all of Bmg can be precomputed, but that is for later
  MATRIX H = computeH(gaussPoint);
  MATRIX DmHg = _Dm * H;
  Real dXdXi = DmHg.determinant();
  MATRIX DmHgInverse = DmHg.inverse();
  MATRIX Bmg = DmHgInverse.transpose() * H.transpose() * dXdXi * weight;

  MATRIX F = computeF(gaussPoint);

  /*
  //MATRIX S = secondPiolaKirchhoffStVK(gaussPoint);
  MATRIX S = _material->PK2(F);

  //MATRIX IS = diag(S);
  MATRIX IS = DFSDF(S);

  MATRIX diagF = blockDiag(F, 2);

  // this can be precomputed once per material setting, at least for StVK
  //MATRIX DSDE = DSDEStVK();
  MATRIX DSDE = _material->DSDE();

  // no getting around computing this one -- changes with F
  MATRIX DEDF = computeDEDF(F);

  MATRIX A = IS + diagF * DSDE * DEDF;
  */
  MATRIX DPDF = _material->DPDF(F);

  // do it per entry
  MATRIX K(8,8);
  for (int i = 0; i < 8; i++)
  {
    MATRIX DFDu_i = computeDFDu(gaussPoint, i);
    VECTOR flatDFDu_i = flatten(DFDu_i);
    VECTOR flatDPDu = DPDF * flatDFDu_i;

    MATRIX DPDu = unflatten(flatDPDu, 2);
    MATRIX squareKi = DPDu * Bmg;

    K.col(i) = flatten(squareKi);
  }

  return -1.0 * K;
}

///////////////////////////////////////////////////////////////////////
// get the forces Jacobian at a Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeForceJacobianFast(const int gaussPoint)
{
  TIMER functionTimer(__FUNCTION__);
  // all of Bmg can be precomputed, but that is for later
  MATRIX& H = _H[gaussPoint];
  MATRIX& Bmg = _Bmg[gaussPoint];

  MATRIX Ds = vertexMatrix();
  MATRIX F = (Ds * H) * _DmH_Inv[gaussPoint];

  /*
  //MATRIX S = secondPiolaKirchhoffStVK(gaussPoint);
  MATRIX S = _material->PK2(F);

  //MATRIX IS = diag(S);
  MATRIX IS = DFSDF(S);

  MATRIX diagF = blockDiag(F, 2);

  // this can be precomputed once per material setting, at least for StVK
  //MATRIX DSDE = DSDEStVK();
  MATRIX DSDE = _material->DSDE();

  // no getting around computing this one -- changes with F
  MATRIX DEDF = computeDEDF(F);

  //MATRIX A = IS + diagF * DSDE * DEDF;
  MATRIX DPDF = IS + diagF * DSDE * DEDF;
  */
  MATRIX DPDF = _material->DPDF(F);
  //cout << " DPDF: " << endl;
  //cout << DPDF << endl;
  //cout << " Dims of DPDF: " << DPDF.rows() << " " << DPDF.cols() << endl;

#if 0
  // indefiniteness fix
  EigenSolver<MATRIX> eigensolver(DPDF);
  MATRIX vectors = eigensolver.eigenvectors().real();

  MATRIX lambda(4,4);
  lambda.setZero();
  for (int i = 0; i < 4; i++)
  {
    Real value = eigensolver.eigenvalues().real()[i];
    lambda(i,i) = (value < 0) ? 0 : value;
  }
  //cout << " eigenvalues: " << eigensolver.eigenvalues().real().transpose() << endl;

  //cout << " before: " << endl << DPDF << endl;
  DPDF = vectors * lambda * vectors.inverse();
  //cout << " after: " << endl << DPDF << endl;
#endif

  // do it per entry
  // it actually spends all its time here, and yet there is no material call
  // on the inside, suggesting that something highly optimized should be
  // possible instead
  TIMER finalBuild("K final build");
  MATRIX K(8,8);
  for (int i = 0; i < 8; i++)
  {
    /*
    MATRIX DFDu_i = computeDFDu(gaussPoint, i);
    //cout << " Dims of DFDu: " << DFDu_i.rows() << " " << DFDu_i.cols() << endl;
    VECTOR flatDFDu_i = flatten(DFDu_i);
    //VECTOR flatDPDu = A * flatDFDu_i;
    VECTOR flatDPDu = DPDF * flatDFDu_i;
    */
    VECTOR flatDPDu = DPDF * _flatDFDu[gaussPoint][i];
    MATRIX DPDu = unflatten(flatDPDu, 2);

    //cout << " Dims of DPDu: " << DPDu.rows() << " " << DPDu.cols() << endl;
    //exit(0);
    MATRIX squareKi = DPDu * Bmg;

    K.col(i) = flatten(squareKi);
  }
  finalBuild.stop();

  // probe the eigenvalues
  // BAD
  //SelfAdjointEigenSolver<MATRIX> eigensolver(DPDF);
  //cout << " eigenvalues: " << eigensolver.eigenvalues().transpose() << endl;

#if 1
  /*
  // indefiniteness fix
  //SelfAdjointEigenSolver<MATRIX> eigensolver(K);
  EigenSolver<MATRIX> eigensolver(K);
  MATRIX vectors = eigensolver.eigenvectors().real();

  MATRIX lambda(8,8);
  lambda.setZero();
  for (int i = 0; i < 8; i++)
  {
    Real value = eigensolver.eigenvalues().real()[i];
    //lambda(i,i) = (value > 0) ? -1e-4 : value;
    lambda(i,i) = (value < 0) ? 0 : value;
  }
  cout << " eigenvalues: " << eigensolver.eigenvalues().real().transpose() << endl;

  cout << " before: " << endl << K << endl;
  K = vectors * lambda * vectors.inverse();
  cout << " after: " << endl << K << endl;
  */
  /*
  // cache the condition number
  // BAD
  SelfAdjointEigenSolver<MATRIX> eigensolver(K);
  MATRIX vectors = eigensolver.eigenvectors();

  _minEig = (eigensolver.eigenvalues()[0]);
  _maxEig = (eigensolver.eigenvalues()[0]);
  for (int i = 1; i < 8; i++)
  {
    if ((eigensolver.eigenvalues()[i]) < _minEig)
      _minEig = (eigensolver.eigenvalues()[i]);
    if (eigensolver.eigenvalues()[i] > _maxEig)
      _maxEig = (eigensolver.eigenvalues()[i]);
  }
  */

  /*
  _minEig = fabs(eigensolver.eigenvalues()[0]);
  _maxEig = fabs(eigensolver.eigenvalues()[0]);
  for (int i = 1; i < 8; i++)
  {
    if (fabs(eigensolver.eigenvalues()[i]) < _minEig)
      _minEig = fabs(eigensolver.eigenvalues()[i]);
    if (fabs(eigensolver.eigenvalues()[i]) > _maxEig)
      _maxEig = fabs(eigensolver.eigenvalues()[i]);
  }

  if (_minEig < 1e-8)
    _minEig = 1e-8;
  _conditionNumber = fabs(_maxEig / _minEig);
  */
#endif

  return -1.0 * K;
}

///////////////////////////////////////////////////////////////////////
// get the forces Jacobian at a Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeInvertibleForceJacobian(const int gaussPoint)
{
  TIMER functionTimer(__FUNCTION__);

  // all of Bmg can be precomputed, but that is for later
  MATRIX& H = _H[gaussPoint];
  MATRIX& Bmg = _Bmg[gaussPoint];

  MATRIX Ds = vertexMatrix();
  MATRIX F = (Ds * H) * _DmH_Inv[gaussPoint];

  // get the invertible material. The pointer IS invertible, right?
  // this will throw an ugly exception if that is not true.
  INVERTIBLE* invertible = dynamic_cast<INVERTIBLE*>(_material);

  // use inversion consistency
  MATRIX DPDF = invertible->DPDF(F, _oldDirection[gaussPoint], _oldClamped[gaussPoint],
                                    _newDirection[gaussPoint], _newClamped[gaussPoint]);

  // do it per entry
  // it actually spends all its time here, and yet there is no material call
  // on the inside, suggesting that something highly optimized should be
  // possible instead
  TIMER finalBuild("K final build");
  MATRIX K(8,8);
  for (int i = 0; i < 8; i++)
  {
    VECTOR flatDPDu = DPDF * _flatDFDu[gaussPoint][i];
    MATRIX DPDu = unflatten(flatDPDu, 2);
    MATRIX squareKi = DPDu * Bmg;

    K.col(i) = flatten(squareKi);
  }
  finalBuild.stop();

  return -1.0 * K;
}

///////////////////////////////////////////////////////////////////////
// get the PK1 Jacobian at a Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeDPDu(const VEC2& gaussPoint, const int column)
{
  // the Gauss weight is one for the 2x2 scheme
  Real weight = 1.0;

  // all of Bmg can be precomputed, but that is for later
  MATRIX H = computeH(gaussPoint);
  MATRIX DmHg = _Dm * H;
  Real dXdXi = DmHg.determinant();
  MATRIX DmHgInverse = DmHg.inverse();
  MATRIX Bmg = DmHgInverse.transpose() * H.transpose() * dXdXi * weight;

  // this was computed in second PK, could be cached
  MATRIX F = computeF(gaussPoint);
  /*
  MATRIX S = _material->PK2(F);

  //MATRIX S = secondPiolaKirchhoffStVK(gaussPoint);

  //MATRIX IS = diag(S);
  MATRIX IS = DFSDF(S);

  MATRIX diagF = blockDiag(F, 2);

  // this can be precomputed once per material setting, at least for StVK
  //MATRIX DSDE = DSDEStVK();
  MATRIX DSDE = _material->DSDE();

  // no getting around computing this one -- changes with F
  MATRIX DEDF = computeDEDF(F);

  MATRIX A = IS + diagF * DSDE * DEDF;
  */
  MATRIX DPDF = _material->DPDF(F);

  // do it per entry
  //MATRIX K(8,8);
  //for (int i = 0; i < 8; i++)
  MATRIX DFDu_i = computeDFDu(gaussPoint, column);
  VECTOR flatDFDu_i = flatten(DFDu_i);

  VECTOR flatDPDu = DPDF * flatDFDu_i;

  MATRIX DPDu = unflatten(flatDPDu, 2);

  return DPDu;
}

///////////////////////////////////////////////////////////////////////
// get the PK2 Jacobian at a Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeDSDu(const VEC2& gaussPoint, const int column)
{
  // this can be precomputed once per material setting, at least for StVK
  //MATRIX DSDE = DSDEStVK();
  MATRIX DSDE = _material->DSDE();

  // this was computed in second PK, could be cached
  MATRIX F = computeF(gaussPoint);

  // no getting around computing this one -- changes with F
  MATRIX DEDF = computeDEDF(F);

  MATRIX A = DSDE * DEDF;

  // do it per entry
  MATRIX DFDu_i = computeDFDu(gaussPoint, column);
  VECTOR flatDFDu_i = flatten(DFDu_i);
  VECTOR flatDPDu = A * flatDFDu_i;

  MATRIX DSDu = unflatten(flatDPDu, 2);

  return DSDu;
}

///////////////////////////////////////////////////////////////////////
// get the first Piola Kirchhoff for a single Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::firstPiolaKirchhoff(const VEC2& gaussPoint)
{
  MATRIX2 F = computeF(gaussPoint);
  return _material->PK1(F);
}

///////////////////////////////////////////////////////////////////////
// get the second Piola Kirchhoff for a single Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::secondPiolaKirchhoff(const VEC2& gaussPoint)
{
  MATRIX2 F = computeF(gaussPoint);
  return _material->PK2(F);
}

/*
///////////////////////////////////////////////////////////////////////
// get the first Piola Kirchhoff for Co-rotational
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::firstPiolaKirchhoffCorotational(const VEC2& gaussPoint)
{
  MATRIX2 F = computeF(gaussPoint);
  //JacobiSVD<MATRIX2> svd(F, ComputeThinU | ComputeThinV);
  JacobiSVD<MATRIX2> svd(F, ComputeFullU | ComputeFullV);

  MATRIX2 U = svd.matrixU();
  MATRIX2 V = svd.matrixV();
  MATRIX2 Sigma;
  Sigma.setZero();
  Sigma.diagonal() = svd.singularValues();

  MATRIX2 S = V * Sigma * V.transpose();
  MATRIX2 R = U * V.transpose();

  // settings from the 2008 cubature paper
  const Real lambda = 1000;
  const Real mu = 5000;
  const MATRIX2 I = MATRIX2::Identity();
  const MATRIX2 SminusI = S - I;
  return R * (2.0 * mu * SminusI  + lambda * SminusI.trace() * I);
}
*/

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR SQUARE::forceVector() const
{
  VECTOR final(8);

  int i = 0;
  for (int col = 0; col < 4; col++)
    for (int row = 0; row < 2; row++, i++)
      final[i] = _forces(row, col);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR SQUARE::getDisplacement() const
{
  VECTOR final(8);

  int i = 0;
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 2; y++, i++)
      final[i] = (*_vertices[x])[y] - _restPose[x][y];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::setDisplacement(VECTOR& u)
{
  assert(u.size() == 8);
  int i = 0;
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 2; y++, i++)
      (*_vertices[x])[y] = _restPose[x][y] + u[i];
}

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR SQUARE::flatten(const MATRIX& A)
{
  VECTOR final(A.rows() * A.cols());

  int i = 0;
  for (int y = 0; y < A.cols(); y++)
    for (int x = 0; x < A.rows(); x++, i++)
      final[i] = A(x,y);

  return final;
}

///////////////////////////////////////////////////////////////////////
// unflatten a vector into a matrix, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::unflatten(const VECTOR& v, const int rows)
{
  // make sure the vector size doesn't leave leftovers
  assert(v.size() % rows == 0);

  const int cols = v.size() / rows;
  MATRIX final(rows, cols);

  int i = 0;
  for (int y = 0; y < cols; y++)
    for (int x = 0; x < rows; x++, i++)
      final(x,y) = v[i]; 

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of deformation gradient F
// with respect to displacement u
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeDFDu(const VEC2& Xi)
{
  MATRIX H = computeH(Xi); 

  // everything except the u-varying component;
  MATRIX A = H * (_Dm * H).inverse();

  // probe each entry with a 1 to build each column
  const int rows = 2;
  const int cols = 4;
  MATRIX delta(rows, cols);
  delta.setZero();

  MATRIX final(4,8);
  final.setZero();
  int i = 0;
  for (int y = 0; y < cols; y++)
    for (int x = 0; x < rows; x++, i++)
    {
      delta(x,y) = 1;

      MATRIX probe = delta * A;
      VECTOR flat = flatten(probe);
      final.col(i) = flat;

      // revert the delta function back
      delta(x,y) = 0;
    }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of deformation gradient F
// with respect to displacement u_i
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeDFDu(const VEC2& Xi, const int index)
{
  MATRIX H = computeH(Xi); 

  // everything except the u-varying component;
  MATRIX A = H * (_Dm * H).inverse();

  // probe each entry with a 1 to build each column
  const int rows = 2;
  const int cols = 4;
  MATRIX delta(rows, cols);
  delta.setZero();

  int col = index / 2;
  int row = index % 2;

  delta(row, col) = 1.0;

  return delta * A;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of deformation gradient F
// with respect to displacement u_i
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeDFDu(const int gaussPoint, const int index)
{
  // everything except the u-varying component;
  MATRIX& A = _DFDuA[gaussPoint];

  // probe each entry with a 1 to build each column
  const int rows = 2;
  const int cols = 4;
  MATRIX delta(rows, cols);
  delta.setZero();

  int col = index / 2;
  int row = index % 2;

  delta(row, col) = 1.0;

  return delta * A;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of E = 1/2 (F^T * F - I)
// with respect to F
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::computeDEDF(const MATRIX2& F)
{
  const Real& f00 = F(0,0);
  const Real& f10 = F(1,0);
  const Real& f01 = F(0,1);
  const Real& f11 = F(1,1);
  MATRIX final(4,4);

  final(0,0) = 2.0 * f00;
  final(1,0) = f01;
  final(2,0) = f01;
  final(3,0) = 0.0;

  final(0,1) = 2.0 * f10;
  final(1,1) = f11;
  final(2,1) = f11;
  final(3,1) = 0.0;

  final(0,2) = 0.0;
  final(1,2) = f00;
  final(2,2) = f00;
  final(3,2) = 2.0 * f01;
  
  final(0,3) = 0.0;
  final(1,3) = f10;
  final(2,3) = f10;
  final(3,3) = 2.0 * f11;

  // E = 0.5 (F^T * F - I), so a 0.5 factor is needed in front
  final *= 0.5;

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of S, the second Piola-Kirchhoff
// with respect to E = 1/2 (F^T * F - I)
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::DSDE()
{
  // this is specifically for StVK
  /*
  MATRIX final(4,4);
  final.setZero();

  // settings from the 2008 cubature paper
  const Real lambda = 1000;
  const Real mu = 5000;

  final(0,0) = lambda + 2.0 * mu;
  final(3,0) = lambda;

  final(1,1) = 2.0 * mu;

  final(2,2) = 2.0 * mu;

  final(0,3) = lambda;
  final(3,3) = lambda + 2.0 * mu;

  return final;
  */
  return _material->DSDE();
}

///////////////////////////////////////////////////////////////////////
// take the gradient of FS (i.e. PK1) w.r.t F assuming that S is frozen
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::DFSDF(const MATRIX& S)
{
  MATRIX final(4,4);
  final.setZero();

  final(0,0) = S(0,0);
  final(1,1) = S(0,0);
  final(2,2) = S(1,1);
  final(3,3) = S(1,1);

  final(2,0) = S(0,1);
  final(3,1) = S(0,1);
  
  final(0,2) = S(1,0);
  final(1,3) = S(1,0);

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the entries of the matrix and put them along the diagonal
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::diag(const MATRIX& A)
{
  const int entries = A.rows() * A.cols();
  MATRIX final(entries, entries);
  final.setZero();

  int i = 0;
  for (int y = 0; y < A.cols(); y++)
    for (int x = 0; x < A.rows(); x++, i++)
      final(i,i) = A(x,y);

  return final;
}

///////////////////////////////////////////////////////////////////////
// repeat the given matrix "repeat" times along the block diagonal
///////////////////////////////////////////////////////////////////////
MATRIX SQUARE::blockDiag(const MATRIX& A, const int repeats)
{
  const int rows = A.rows();
  const int cols = A.cols();
  MATRIX final(rows * repeats, cols * repeats);
  final.setZero();

  for (int y = 0; y < cols; y++)
    for (int x = 0; x < rows; x++)
    {
      for (int i = 0; i < repeats; i++)
        final(i * rows + x, i * cols + y) = A(x,y);
    }

  return final;
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::squash()
{
#if 1
  //scale(VEC2(2, 1));
  scale(VEC2(3, 5));
#else
 /* 
  scale(VEC2(1.1, 1));
  scale(VEC2(1.3, 1));
  //scale(VEC2(0.5, 1.5));
  scale(VEC2(3.14159, 2.65));
  rotate(M_PI / 8);
  (*_vertices[0])[0] += 0.1;
  */
#endif
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::testK()
{
  VECTOR uOriginal = getDisplacement();

  Real eps = 1e-2;

  MATRIX K = computeForceJacobian();

  cout << "======================================" << endl;
  cout << " Testing K (DfDu): " << endl;
  cout << "======================================" << endl;

  Real diff = 0;
  MATRIX Knumerical(8,8);
  for (int i = 0; i < 6; i++)
  {
    for (int x = 0; x < 8; x++)
    {
      setDisplacement(uOriginal);
      computeForces();
      VECTOR force0 = forceVector();

      VECTOR uDelta = uOriginal;
      uDelta[x] += eps;
      setDisplacement(uDelta);
      computeForces();
      VECTOR force1 = forceVector();

      VECTOR column = force1 - force0;
      column *= 1.0 / eps;

      Knumerical.col(x) = column;
    }

    MATRIX Kdiff = K - Knumerical;
    diff = Kdiff.norm() / Kdiff.size();
    cout << " eps: \t" << eps << " \t diff: \t" << diff << endl;

    eps *= 0.1;
  }
  if (diff < 1e-4)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::testRestStability()
{
  MATRIX2 I;
  I.setZero();
  I(0,0) = I(1,1) = 1.0;

  MATRIX P = _material->PK1(I);
  Real norm = P.squaredNorm();

  cout << "======================================" << endl;
  cout << " Testing rest stability " << endl;
  cout << "======================================" << endl;
  cout << " Squared norm: " << norm << endl;

  cout << "P: " << endl;
  cout << P << endl;

  if (norm < 1e-7)
    cout << " TEST PASSED " << endl;
  else
  {
    cout << " TEST FAILED " << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::testDPDu()
{
  VECTOR uOriginal = getDisplacement();

  //Real eps = 1e-7;
  Real eps = 1e-2;

  const Real invSqrt3 = 1.0 / sqrt(3.0);
  VEC2 Xi(1, 1);
  Xi *= invSqrt3;

  MATRIX DPDu(4,8);
  for (int i = 0; i < 8; i++)
  {
    MATRIX col = computeDPDu(Xi, i);
    DPDu.col(i) = SQUARE::flatten(col);
  }
  setDisplacement(uOriginal);
  MATRIX PK10 = firstPiolaKirchhoff(Xi);

  cout << "======================================" << endl;
  cout << " Testing DPDu: " << endl;
  cout << "======================================" << endl;

  Real diff = 0;
  MATRIX Pnumerical(4,8);
  for (int i = 0; i < 6; i++)
  {
    for (int x = 0; x < 8; x++)
    {
      VECTOR uDelta = uOriginal;
      uDelta[x] += eps;
      setDisplacement(uDelta);
      MATRIX PK11 = firstPiolaKirchhoff(Xi);

      MATRIX diff = PK11 - PK10;
      diff *= 1.0 / eps;

      Pnumerical.col(x) = SQUARE::flatten(diff);
    }

    MATRIX Pdiff = DPDu - Pnumerical;
    diff = Pdiff.norm() / Pdiff.size();
    cout << " eps: \t" << eps << " \t diff: \t" << diff << endl;

#if 1
    if (i == 5)
    {
      cout << " direct: " << endl;
      cout << DPDu << endl;
      cout << " numerical: " << endl;
      cout << Pnumerical << endl;
      cout << " diff: " << endl;
      cout << Pdiff << endl;
    }
#endif

    eps *= 0.1;
  }
  if (diff < 1e-4)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::testDpsiDF()
{
  VECTOR uOriginal = getDisplacement();

  //Real eps = 1e-7;
  Real eps = 1e-2;

  const Real invSqrt3 = 1.0 / sqrt(3.0);
  VEC2 Xi(1, 1);
  Xi *= invSqrt3;
  MATRIX F = computeF(Xi);

  MATRIX DpsiDF = _material->PK1(F);
  Real psi = _material->psi(F);

  cout << "Psi: " << psi << endl;
  cout << "P: " << endl;
  cout << DpsiDF << endl;
  cout << "F: " << F << endl;

  setDisplacement(uOriginal);

  cout << "======================================" << endl;
  cout << " Testing DpsiDF: " << endl;
  cout << "======================================" << endl;

  Real diff = 0;
  MATRIX Pnumerical(2,2);
  for (int e = 0; e < 6; e++)
  {
    MATRIX2 F0 = F;
    Real P0 = _material->psi(F0);

    int i = 0;
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 2; x++, i++)
      {
        MATRIX2 F1 = F;
        F1(x,y) += eps;

        Real P1 = _material->psi(F1);

        Real diff = P1 - P0;
        diff *= 1.0 / eps;

        Pnumerical(x,y) = diff;
      }

    MATRIX Pdiff = DpsiDF - Pnumerical;
    diff = Pdiff.norm() / Pdiff.size();

#if 1
    if (e == 5)
    {
      cout << " numerical: " << endl;
      cout << Pnumerical << endl;
      cout << " diff: " << endl;
      cout << Pdiff << endl;
    }
#endif
    cout << " eps: \t" << eps << " \t diff: \t" << diff << endl;

    eps *= 0.1;
  }

  if (diff < 1e-3)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::testDPDF()
{
  VECTOR uOriginal = getDisplacement();

  //Real eps = 1e-7;
  Real eps = 1e-2;

  const Real invSqrt3 = 1.0 / sqrt(3.0);
  VEC2 Xi(1, 1);
  Xi *= invSqrt3;
#if 1
  MATRIX F = computeF(Xi);
#else
  MATRIX F(2,2);
  F.setZero();
  //F(0,0) = 1.5;
  F(0,0) = 1.0;
  F(1,1) = 1.0;
#endif

  MATRIX DPDF = _material->DPDF(F);

  cout << "DPDF: " << endl;
  cout << DPDF << endl;
  cout << "Psi: " << endl;
  cout << _material->psi(F) << endl;
  cout << "F: " << F << endl;

  /*
  cout << " Mu: " << ((HYPER_TAN*)_material)->mu(F) << endl;
  cout << " P: " << endl;
  cout << _material->PK1(F) << endl;
  */

  setDisplacement(uOriginal);
  //MATRIX PK10 = firstPiolaKirchhoff(Xi);
  //MATRIX PK10 = firstPiolaKirchhoff(Xi);

  cout << "======================================" << endl;
  cout << " Testing DPDF: " << endl;
  cout << "======================================" << endl;

  Real diff = 0;
  MATRIX Pnumerical(4,4);
  for (int e = 0; e < 6; e++)
  {
    MATRIX2 F0 = F;
    MATRIX2 P0 = _material->PK1(F0);

    int i = 0;
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 2; x++, i++)
      {
        MATRIX2 F1 = F;
        F1(x,y) += eps;

        MATRIX2 P1 = _material->PK1(F1);

        MATRIX diff = P1 - P0;
        diff *= 1.0 / eps;

        Pnumerical.col(i) = SQUARE::flatten(diff);
      }

    MATRIX Pdiff = DPDF - Pnumerical;
    diff = Pdiff.norm() / Pdiff.size();

#if 1
    if (e == 5)
    {
      cout << " numerical: " << endl;
      cout << Pnumerical << endl;
      cout << " diff: " << endl;
      cout << Pdiff << endl;
    }
#endif
    cout << " eps: \t" << eps << " \t diff: \t" << diff << endl;

    eps *= 0.1;
  }

  if (diff < 1e-3)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::testDRDF()
{
  //VECTOR uOriginal = getDisplacement();

  //Real eps = 1e-7;
  Real eps = 1e-2;

  const Real invSqrt3 = 1.0 / sqrt(3.0);
  VEC2 Xi(1, 1);
  Xi *= invSqrt3;
  MATRIX F = computeF(Xi);

  // this test is only valid for corotational, so if this cast fails,
  // you're doing it wrong.
  COROTATIONAL* corotational = (COROTATIONAL*)_material;
  MATRIX DRDF = corotational->DRDF(F);
  DRDF *= -1.0;

  cout << "DRDF: " << endl;
  cout << DRDF << endl;
  //exit(0);

  //setDisplacement(uOriginal);
  MATRIX2 R0 = corotational->R(F);

  cout << "======================================" << endl;
  cout << " Testing DRDF: " << endl;
  cout << "======================================" << endl;

  cout << " Initial S: " << endl;
  cout << corotational->S(F) << endl;

  Real diff = 0;
  MATRIX Rnumerical(4,4);
  for (int e = 0; e < 6; e++)
  {
    int i = 0;
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 2; x++, i++)
      {
        MATRIX2 F1 = F;
        F1(x,y) += eps;

        MATRIX2 R1 = corotational->R(F1);

        MATRIX2 diff = R1 - R0;
        diff *= 1.0 / eps;

        Rnumerical.col(i) = SQUARE::flatten(diff);
      }

    MATRIX Rdiff = DRDF - Rnumerical;
    diff = Rdiff.norm() / Rdiff.size();

    cout << " eps: \t" << eps << " \t diff: \t" << diff << endl;
#if 1
    if (e == 5)
    {
      cout << " numerical: " << endl;
      cout << Rnumerical << endl;
      cout << " diff: " << endl;
      cout << Rdiff << endl;
    }
#endif

    eps *= 0.1;
  }

  if (diff < 1e-3)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::testDSDF()
{
  //VECTOR uOriginal = getDisplacement();

  //Real eps = 1e-7;
  Real eps = 1e-2;

  const Real invSqrt3 = 1.0 / sqrt(3.0);
  VEC2 Xi(1, 1);
  Xi *= invSqrt3;
  MATRIX F = computeF(Xi);

  // this test is only valid for corotational, so if this cast fails,
  // you're doing it wrong.
  COROTATIONAL* corotational = (COROTATIONAL*)_material;
  MATRIX DSDF = corotational->DSDF(F);

  cout << "DSDF: " << endl;
  cout << DSDF << endl;
  //exit(0);

  //setDisplacement(uOriginal);
  MATRIX2 S0 = corotational->S(F);

  cout << "======================================" << endl;
  cout << " Testing DSDF: " << endl;
  cout << "======================================" << endl;

  Real diff = 0;
  MATRIX Snumerical(4,4);
  Snumerical.setZero();
  for (int e = 0; e < 6; e++)
  {
    int i = 0;
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 2; x++, i++)
      {
        MATRIX2 F1 = F;
        F1(x,y) += eps;

        MATRIX2 S1 = corotational->S(F1);
        //cout << " S1: " << endl;
        //cout << S1 << endl;
        //cout << " S0: " << endl;
        //cout << S0 << endl;

        MATRIX2 matrixDiff = S1 - S0;
        //cout << " matrix diff: " << endl;
        //cout << matrixDiff << endl;

        matrixDiff *= 1.0 / eps;

        Snumerical.col(i) = SQUARE::flatten(matrixDiff);
      }

    MATRIX Sdiff = DSDF - Snumerical;
    diff = Sdiff.norm() / Sdiff.size();

    if (e == 5)
    {
      cout << " numerical: " << endl;
      cout << Snumerical << endl;
      cout << " diff: " << endl;
      cout << Sdiff << endl;
    }
    cout << " eps: \t" << eps << " \t diff: \t" << diff << endl;

    eps *= 0.1;
  }

  if (diff < 1e-3)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::testDSDu()
{
  VECTOR uOriginal = getDisplacement();

  //Real eps = 1e-7;
  Real eps = 1e-2;

  const Real invSqrt3 = 1.0 / sqrt(3.0);
  VEC2 Xi(1, 1);
  Xi *= invSqrt3;

  MATRIX DSDu(4,8);
  for (int i = 0; i < 8; i++)
  {
    MATRIX col = computeDSDu(Xi, i);
    DSDu.col(i) = SQUARE::flatten(col);
  }
  setDisplacement(uOriginal);
  MATRIX PK20 = secondPiolaKirchhoff(Xi);

  cout << "======================================" << endl;
  cout << " Testing DSDu: " << endl;
  cout << "======================================" << endl;

  Real diff = 0;
  MATRIX Snumerical(4,8);
  for (int i = 0; i < 6; i++)
  {
    for (int x = 0; x < 8; x++)
    {
      VECTOR uDelta = uOriginal;
      uDelta[x] += eps;
      setDisplacement(uDelta);
      MATRIX PK21 = secondPiolaKirchhoff(Xi);

      MATRIX diff = PK21 - PK20;
      diff *= 1.0 / eps;

      Snumerical.col(x) = SQUARE::flatten(diff);
    }

    MATRIX Sdiff = DSDu - Snumerical;
    diff = Sdiff.norm() / Sdiff.size();
    cout << " eps: \t" << eps << " \t diff: \t" << diff << endl;

    eps *= 0.1;
  }
  if (diff < 1e-4)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::testDFDu()
{
  VECTOR uOriginal = getDisplacement();

  Real eps = 1e-2;

  const Real invSqrt3 = 1.0 / sqrt(3.0);
  VEC2 Xi(1, 1);
  Xi *= invSqrt3;

  MATRIX DFDu = computeDFDu(Xi);

  cout << "======================================" << endl;
  cout << " Testing DFDu: " << endl;
  cout << "======================================" << endl;

  Real diff = 0;
  MATRIX Fnumerical(4,8);
  for (int i = 0; i < 6; i++)
  {
    for (int x = 0; x < 8; x++)
    {
      setDisplacement(uOriginal);
      MATRIX F0 = computeF(Xi);

      VECTOR uDelta = uOriginal;
      uDelta[x] += eps;
      setDisplacement(uDelta);
      MATRIX F1 = computeF(Xi);

      MATRIX diff = F1 - F0;
      diff *= 1.0 / eps;

      Fnumerical.col(x) = SQUARE::flatten(diff);
    }

    MATRIX Fdiff = DFDu - Fnumerical;
    diff = Fdiff.norm() / Fdiff.size();
    cout << " eps: \t" << eps << " \t diff: \t" << diff << endl;

    eps *= 0.1;
  }
  if (diff < 1e-4)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::testDEDF()
{
  VECTOR uOriginal = getDisplacement();

  Real eps = 1e-2;

  const Real invSqrt3 = 1.0 / sqrt(3.0);
  VEC2 Xi(1, 1);
  Xi *= invSqrt3;

  MATRIX2 F = computeF(Xi);

  MATRIX DEDF = computeDEDF(F);

  cout << "======================================" << endl;
  cout << " Testing DEDF: " << endl;
  cout << "======================================" << endl;

  Real diff = 0;
  MATRIX Enumerical(4,4);
  for (int e = 0; e < 6; e++)
  {
    MATRIX2 F0 = F;
    MATRIX2 E0 = 0.5 * (F0.transpose() * F0 - MATRIX2::Identity());

    int i = 0;
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 2; x++, i++)
      {
        MATRIX2 F1 = F;
        F1(x,y) += eps;

        MATRIX2 E1 = 0.5 * (F1.transpose() * F1 - MATRIX2::Identity());

        MATRIX diff = E1 - E0;
        diff *= 1.0 / eps;

        Enumerical.col(i) = SQUARE::flatten(diff);
      }

    MATRIX Ediff = DEDF - Enumerical;
    diff = Ediff.norm() / Ediff.size();
    cout << " eps: \t" << eps << " \t diff: \t" << diff << endl;

    eps *= 0.1;
  }
  if (diff < 1e-4)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    exit(0);
  }
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::testDEDu()
{
  VECTOR uOriginal = getDisplacement();

  Real eps = 1e-2;

  const Real invSqrt3 = 1.0 / sqrt(3.0);
  VEC2 Xi(1, 1);
  Xi *= invSqrt3;

  MATRIX2 F = computeF(Xi);

  MATRIX DEDF = computeDEDF(F);
  MATRIX DFDu = computeDFDu(Xi);
  MATRIX DEDu = DEDF * DFDu;

  cout << "======================================" << endl;
  cout << " Testing DEDu: " << endl;
  cout << "======================================" << endl;
  //cout << " DEDu: " << endl;
  //cout << DEDu << endl;

  Real diff = 0;
  MATRIX Enumerical(4,8);
  for (int e = 0; e < 6; e++)
  {
    MATRIX2 F0 = F;
    MATRIX2 E0 = 0.5 * (F0.transpose() * F0 - MATRIX2::Identity());

    for (int x = 0; x < 8; x++)
    {
      VECTOR u1 = uOriginal;
      u1[x] += eps;
      setDisplacement(u1);

      MATRIX F1 = computeF(Xi);
      MATRIX2 E1 = 0.5 * (F1.transpose() * F1 - MATRIX2::Identity());

      MATRIX diff = E1 - E0;
      diff *= 1.0 / eps;

      Enumerical.col(x) = SQUARE::flatten(diff);
    }

    MATRIX Ediff = DEDu - Enumerical;
    diff = Ediff.norm() / Ediff.size();
    cout << " eps: \t" << eps << " \t diff: \t" << diff << endl;

    eps *= 0.1;
  }
  if (diff < 1e-4)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////
// run all the unit tests on a deformed square
///////////////////////////////////////////////////////////////////////
void SQUARE::runHypertanUnitTests(MATERIAL* material)
{
  vector<VEC2> vertices(4);
  vertices[0] = VEC2(-1, -1);
  vertices[1] = VEC2(1, -1);
  vertices[2] = VEC2(1, 1);
  vertices[3] = VEC2(-1, 1);

  vector<VEC2*> vertexPointers(4);
  vertexPointers[0] = &vertices[0];
  vertexPointers[1] = &vertices[1];
  vertexPointers[2] = &vertices[2];
  vertexPointers[3] = &vertices[3];

  SQUARE square(material, vertexPointers);
  square.squash();
  //square.stretch();

  //square.testDpsiDF();
  square.testMu();
  square.testDPDF();
  square.testDPDu();
  square.testK();
}

///////////////////////////////////////////////////////////////////////
// run all the unit tests on a deformed square
///////////////////////////////////////////////////////////////////////
void SQUARE::runUnitTests(MATERIAL* material)
{
  vector<VEC2> vertices(4);
  vertices[0] = VEC2(-1, -1);
  vertices[1] = VEC2(1, -1);
  vertices[2] = VEC2(1, 1);
  vertices[3] = VEC2(-1, 1);

  vector<VEC2*> vertexPointers(4);
  vertexPointers[0] = &vertices[0];
  vertexPointers[1] = &vertices[1];
  vertexPointers[2] = &vertices[2];
  vertexPointers[3] = &vertices[3];

  SQUARE square(material, vertexPointers);
  square.squash();
  cout << " ======================================== " << endl;
  cout << " Testing " << material->name().c_str() << " material " << std::endl;
  cout << " ======================================== " << endl;

  square.testDpsiDF();
  
  ///square.testMu();
  /*
  square.testDEDu();
  square.testDEDF();
  square.testDFDu();
  square.testDSDu();
  */
  square.testDPDF();
  square.testDPDu();
  square.testK();
  square.testRestStability();
}

///////////////////////////////////////////////////////////////////////
// run all the unit tests on a deformed square
///////////////////////////////////////////////////////////////////////
void SQUARE::runEigenvalueTest(MATERIAL* material)
{
  vector<VEC2> vertices(4);
  vertices[0] = VEC2(-1, -1);
  vertices[1] = VEC2(1, -1);
  vertices[2] = VEC2(1, 1);
  vertices[3] = VEC2(-1, 1);

  vector<VEC2*> vertexPointers(4);
  vertexPointers[0] = &vertices[0];
  vertexPointers[1] = &vertices[1];
  vertexPointers[2] = &vertices[2];
  vertexPointers[3] = &vertices[3];

  SQUARE square(material, vertexPointers);
  square.squash();
  cout << " ======================================== " << endl;
  cout << " Testing eigenvalue expressions of " << material->name().c_str() << " material " << std::endl;
  cout << " ======================================== " << endl;

  square.testEigenvalues();
  square.testEigensystem();
}

///////////////////////////////////////////////////////////////////////
// see if eigenvalues match up to the call in "eigenvalues"
///////////////////////////////////////////////////////////////////////
void SQUARE::testEigenvalues()
{
  // get the deformation gradient
  const Real invSqrt3 = 1.0 / sqrt(3.0);
  VEC2 Xi(1, 1);
  Xi *= invSqrt3;
  MATRIX F = computeF(Xi);

  // get the Hessian
  MATRIX DPDF = _material->DPDF(F);

  Eigen::SelfAdjointEigenSolver<MATRIX> eigensolver(DPDF);
  if (eigensolver.info() != Eigen::Success) abort();

  VECTOR numerical = eigensolver.eigenvalues();
  cout << "The numerical eigenvalues are:\n" << numerical << endl;

  VECTOR direct = _material->eigenvalues(F);
  cout << "The direct eigenvalues are:\n" << direct << endl;

  VECTOR diff = numerical - direct;
  Real diffNorm = diff.norm() / numerical.norm();
  cout << " Diff magnitude: " << diffNorm << endl;

  if (diffNorm > 1e-6)
  {
    cout << " TEST FAILED " << endl;
    exit(0);
  }
  cout << " TEST PASSED " << endl;
}

///////////////////////////////////////////////////////////////////////
// see if eigenvalues match up to the call in "eigenvectors"
///////////////////////////////////////////////////////////////////////
void SQUARE::testEigensystem()
{
  // get the deformation gradient
  const Real invSqrt3 = 1.0 / sqrt(3.0);
  VEC2 Xi(1, 1);
  Xi *= invSqrt3;
  MATRIX F = computeF(Xi);

  // get the Hessian
  MATRIX DPDF = _material->DPDF(F);

  Eigen::SelfAdjointEigenSolver<MATRIX> eigensolver(DPDF);
  if (eigensolver.info() != Eigen::Success) abort();

  VECTOR numericalValues = eigensolver.eigenvalues();
  cout << "The numerical eigenvalues are:\n" << endl << numericalValues << endl;

  MATRIX numericalVectors = eigensolver.eigenvectors();
  cout << "The numerical eigenvectors are:\n" << endl << numericalVectors << endl;

  VEC4 directValues;
  MATRIX4 directVectors;
  _material->eigensystem(F, directVectors, directValues);
  cout << "The direct eigenvalues are:\n" << directValues << endl;
  cout << "The direct eigenvectors are:\n" << directVectors << endl;

  // see that the eigenvalues check out
  VECTOR diffValues = numericalValues - directValues;
  Real diffValuesNorm = diffValues.norm() / numericalValues.norm();
  cout << " Diff magnitude: " << diffValuesNorm << endl;
  if (diffValuesNorm > 1e-6)
  {
    cout << " EIGENVALUE TEST FAILED " << endl;
    exit(0);
  }
  cout << " EIGENVALUE TEST PASSED " << endl;

  // see that the eigenvectors check out
  MATRIX4 product = directVectors.transpose() * numericalVectors;
  cout << " The product is:\n " << product << endl;

  // square the product to get rid of the signs
  MATRIX4 diffVectors = product * product - MATRIX4::Identity();
  Real diffVectorsNorm = diffVectors.norm();
  cout << " Diff magnitude: " << diffVectorsNorm << endl;
  if (diffValuesNorm > 1e-6)
  {
    cout << " EIGENVECTOR TEST FAILED " << endl;
    exit(0);
  }
  cout << " EIGENVECTOR TEST PASSED " << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::testMu()
{
  VECTOR uOriginal = getDisplacement();

  //Real eps = 1e-7;
  Real eps = 1e-2;

  const Real invSqrt3 = 1.0 / sqrt(3.0);
  VEC2 Xi(1, 1);
  Xi *= invSqrt3;
  //MATRIX F = computeF(Xi);
  MATRIX F(2,2);
  F.setZero();
  F(0,0) = 1.5;
  F(1,1) = 1.0;

  HYPER_TAN* hypertan = (HYPER_TAN*)_material;

  MATRIX DmuDF = hypertan->DmuDF(F);
  Real mu = hypertan->psi(F);

  cout << "mu: " << mu << endl;
  cout << "P: " << hypertan->PK1(F) << endl;
  cout << "DmuDF: " << endl;
  cout << DmuDF << endl;

  setDisplacement(uOriginal);

  cout << "======================================" << endl;
  cout << " Testing DmuDF: " << endl;
  cout << "======================================" << endl;

  Real diff = 0;
  MATRIX munumerical(2,2);
  for (int e = 0; e < 6; e++)
  {
    MATRIX2 F0 = F;
    Real mu0 = hypertan->mu(F0);

    int i = 0;
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 2; x++, i++)
      {
        MATRIX2 F1 = F;
        F1(x,y) += eps;

        Real mu1 = hypertan->mu(F1);

        Real diff = mu1 - mu0;
        diff *= 1.0 / eps;

        munumerical(x,y) = diff;
      }

    MATRIX mudiff = DmuDF - munumerical;
    diff = mudiff.norm() / mudiff.size();

#if 1
    if (e == 5)
    {
      cout << " numerical: " << endl;
      cout << munumerical << endl;
      cout << " diff: " << endl;
      cout << mudiff << endl;
    }
#endif
    cout << " eps: \t" << eps << " \t diff: \t" << diff << endl;

    eps *= 0.1;
  }

  if (diff < 1e-3)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::testSVD()
{
  const Real invSqrt3 = 1.0 / sqrt(3.0);
  VEC2 Xi(1, 1);
  Xi *= invSqrt3;
  MATRIX F = computeF(Xi);

  VECTOR Svec;
  MATRIX2 U, S, V;
  COROTATIONAL* corotational = (COROTATIONAL*)_material;
  corotational->svd(F, U, Svec, V);

  S.setZero();
  S(0,0) = Svec[0];
  S(1,1) = Svec[1];

  /*
  cout << " U: " << endl;
  cout << U << endl;
  cout << " S: " << endl;
  cout << S << endl;
  cout << " V: " << endl;
  cout << V << endl;
  */

  cout << "======================================" << endl;
  cout << " Testing SVD gradient: " << endl;
  cout << "======================================" << endl;
  Real totalDiff = 0;
  for (int j = 0; j < 2; j++)
    for (int i = 0; i < 2; i++)
    {
      //int i = 0;
      //int j = 1;
      MATRIX2 DU, DS, DV;
      corotational->svdGradient(F, i, j, DU, DS, DV);

      /*
      cout << " DU: " << endl;
      cout << DU << endl;
      cout << " DS: " << endl;
      cout << DS << endl;
      cout << " DV: " << endl;
      cout << DV << endl;
      */

      MATRIX2 product = DU * S * V.transpose() + U * DS * V.transpose() + U * S * DV.transpose();

      MATRIX2 delta;
      delta.setZero();
      delta(i,j) = 1;
      //cout << delta << endl;

      //cout << " product: " << endl;
      //cout << product << endl;
      MATRIX2 gDiff = delta - product;
      Real diff = gDiff.norm();
      cout << " i,j " << i << ", " << j << " diff: " << diff << endl;

      totalDiff += diff;
    }
  if (totalDiff < 1e-8)
    cout << " TEST PASSED " << endl;
  else
  {
    cout << " TEST FAILED " << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////
// run all the unit tests on a deformed square
///////////////////////////////////////////////////////////////////////
void SQUARE::runCorotationalUnitTests(COROTATIONAL* corotational)
{
  vector<VEC2> vertices(4);
  vertices[0] = VEC2(-1, -1);
  vertices[1] = VEC2(1, -1);
  vertices[2] = VEC2(1, 1);
  vertices[3] = VEC2(-1, 1);

  vector<VEC2*> vertexPointers(4);
  vertexPointers[0] = &vertices[0];
  vertexPointers[1] = &vertices[1];
  vertexPointers[2] = &vertices[2];
  vertexPointers[3] = &vertices[3];

  SQUARE square(corotational, vertexPointers);

  /*
  // first run on identity, since the polar gradient can blow up
  // there as well
  cout << "=============================================================== " << endl;
  cout << " Running on identity" << endl;
  cout << "=============================================================== " << endl;
  //square.testDRDF();
  //square.testDSDF();
  square.testDPDF();
  square.testK();
  */

  // then run with it squashed
  cout << "=============================================================== " << endl;
  cout << " Running on squashed " << endl;
  cout << "=============================================================== " << endl;
  square.squash();
  square.testDRDF();
  square.testDSDF();
  exit(0);
  MATRIX2 F = square.gaussF();
  //corotational->testPolarGrad(F);
  square.testDpsiDF();
  //square.testDSDF();
  //square.testSVD();
  square.testDPDF();
  square.testK();
}

///////////////////////////////////////////////////////////////////////
// compute the area of the square by cutting into two triangles
// and the using Heron's formula:
//
// http://mathworld.wolfram.com/HeronsFormula.html
///////////////////////////////////////////////////////////////////////
Real SQUARE::area() const
{
  Real final = 0;

  Real a, b, c, s;
  a = (*_vertices[0] - *_vertices[1]).norm();
  b = (*_vertices[1] - *_vertices[2]).norm();
  c = (*_vertices[2] - *_vertices[0]).norm();
  s = 0.5 * (a + b + c);

  // area of the first triangle
  final += sqrt(s * (s - a) * (s - b) * (s - c));

  // keep c, do the other two over
  a = (*_vertices[0] - *_vertices[3]).norm();
  b = (*_vertices[3] - *_vertices[2]).norm();
  s = 0.5 * (a + b + c);
 
  // area of the second triangle 
  final += sqrt(s * (s - a) * (s - b) * (s - c));

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the strain energy inside the square
///////////////////////////////////////////////////////////////////////
Real SQUARE::psi() const
{
  const Real weight = 1.0;
  const Real invSqrt3 = 1.0 / sqrt(3.0);
  const VEC2 Xi[] = {invSqrt3 * VEC2(-1,-1),
                     invSqrt3 * VEC2(1,-1),
                     invSqrt3 * VEC2(1,1),
                     invSqrt3 * VEC2(-1,1)};

  Real final = 0;
  for (int x = 0; x < 4; x++)
    final += _material->psi(computeF(Xi[x]));

  return final;
}

///////////////////////////////////////////////////////////////////////
// color ramp function (used to show the condition number)
///////////////////////////////////////////////////////////////////////
VEC3 SQUARE::ramp(float value)
{
  const int NUM_COLORS = 4;
  static float color[NUM_COLORS][3] = //{ {0,0,1}, {0,1,0}, {1,1,0}, {1,0,0} };
  {
  {254 / 255.0,240 / 255.0,217 / 255.0},{253 / 255.0,204 / 255.0,138 / 255.0},{252 / 255.0,141 / 255.0,89 / 255.0},{215 / 255.0,48 / 255.0, 31 / 255.0}
  };

  // A static array of 4 colors:  (blue,   green,  yellow,  red) using {r,g,b} for each.
  int idx1;        // |-- Our desired color will be between these two indexes in "color".
  int idx2;        // |
  float fractBetween = 0;  // Fraction between "idx1" and "idx2" where our value is.

  if (value <= 0)      {  idx1 = idx2 = 0;            }    // accounts for an input <=0
  else if(value >= 1)  {  idx1 = idx2 = NUM_COLORS-1; }    // accounts for an input >=0
  else
  {
    value = value * (NUM_COLORS-1);        // Will multiply value by 3.
    idx1  = floor(value);                  // Our desired color will be after this index.
    idx2  = idx1+1;                        // ... and before this index (inclusive).
    fractBetween = value - float(idx1);    // Distance between the two indexes (0-1).
  }

  VEC3 final;

  if (isnan(value))
  {
    final.setZero();
    return final;
  }
  final[0] = (color[idx2][0] - color[idx1][0])*fractBetween + color[idx1][0];
  final[1] = (color[idx2][1] - color[idx1][1])*fractBetween + color[idx1][1];
  final[2] = (color[idx2][2] - color[idx1][2])*fractBetween + color[idx1][2];

  return final;
}

///////////////////////////////////////////////////////////////////////
// run some SVD-related tests
///////////////////////////////////////////////////////////////////////
void SQUARE::svdTest()
{
  // build the Gauss points in the canonical volume
  const Real invSqrt3 = 1.0 / sqrt(3.0);
  const VEC2 Xi[] = {invSqrt3 * VEC2(-1,-1),
                     invSqrt3 * VEC2(1,-1),
                     invSqrt3 * VEC2(1,1),
                     invSqrt3 * VEC2(-1,1)};

  //for (int x = 0; x < 4; x++)
  int x = 3;
  {
    MATRIX H = computeH(Xi[x]);
    MATRIX Ds = vertexMatrix();

    MATRIX2 DmH_Inv = (_Dm * H).inverse();
    MATRIX2 F = (Ds * H) * DmH_Inv;

    cout << " F: " << endl << F << endl;
    //cout << " det F: " << F.determinant() << endl;
    JacobiSVD<MATRIX> svd(F, ComputeFullU | ComputeFullV);
    MATRIX2 U = svd.matrixU();
    MATRIX2 V = svd.matrixV();
    VECTOR Svec = svd.singularValues();

    MATRIX2 S;
    S.setZero();
    S(0,0) = Svec[0];
    S(1,1) = Svec[1];
    cout << " U: " << endl << U << endl;
    cout << " S: " << endl << S << endl;
    cout << " V: " << endl << V << endl;

    MATRIX2 L;
    L.setZero();
    //L(0,0) = 1.0;
    //L(1,1) = (V * U.transpose()).determinant();
    L(0,0) = 1;
    L(1,1) = (U * V.transpose()).determinant();

    //MATRIX2 R = V * L * U.transpose();
    MATRIX2 R = U * L * V.transpose();

    MATRIX2 Unew = U * L;
    MATRIX2 Snew = S * L;
    MATRIX2 Fnew = Unew * Snew * V.transpose();
    cout << " Fnew: " << endl << Fnew << endl;

    //cout << " R: " << endl << U * V.transpose() << endl;
    //cout << " R: " << endl << V * U.transpose() << endl;
    cout << " R: " << endl << R << endl;
  }
}

///////////////////////////////////////////////////////////////////////
// update the inversion direction information
///////////////////////////////////////////////////////////////////////
void SQUARE::updateInversionDirections()
{
  for (int x = 0; x < 4; x++)
  {
    _oldDirection[x] = _newDirection[x];
    _oldClamped[x] = _newClamped[x];

    _drawDirection[x] = _newDirection[x];
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void SQUARE::printReflectionInformation()
{
  cout << "====================================== " << endl;
  cout << " Reflection directions: " << endl;
  cout << "====================================== " << endl;

  for (int x = 0; x < 4; x++)
  {
    cout << "old: " << endl << _oldDirection[x] << endl;
    cout << "new: " << endl << _newDirection[x] << endl;
  }
}
