#include "TRIANGLE.h"
#include <GLUT/glut.h>

// this is needed for a unit test
#include "ARAP.h"

#include <iostream>
using namespace std;

//#define SUPER_GLUE 1000.0
//#define SUPER_GLUE 100.0    // usual setting
//#define SUPER_GLUE 10.0
//#define SUPER_GLUE 1.0

//#define SUPER_GLUE 0.0

#define USING_DS_ZEROING 1
#define USING_ANISOTROPIC_ARAP 0

//Real TRIANGLE::superGlue = 0.0;
Real TRIANGLE::superGlue = 100.0;

TRIANGLE::TRIANGLE(MATERIAL* material, const vector<VEC2*>& vertices) :
  _material(material)
{
  assert(vertices.size() == 3);
  _vertices[0] = vertices[0];
  _vertices[1] = vertices[1];
  _vertices[2] = vertices[2];
  //_vertices[2] = VEC2(1, 1);
  _Dm = vertexMatrix();

  // store these as the rest pose, for now
  for (unsigned int x = 0; x < 3; x++)
    _restPose[x] = *_vertices[x];

  //_forces.resize(2,3);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE::applyF()
{
  MATRIX2 F = computeF();

  for (int x = 0; x < 3; x++)
    (*_vertices[x]) = F * _restPose[x];
}

/*
///////////////////////////////////////////////////////////////////////
// Draw the square to GL
///////////////////////////////////////////////////////////////////////
void TRIANGLE::draw(const VEC3& color)
{
  glLineWidth(5);
  glColor4f(color[0], color[1], color[2], 1.0);
  Real fade = 1;

  glBegin(GL_LINE_STRIP);
    for (int x = 0; x < 3; x++)
    {
      glVertex2f(_vertices[x][0], _vertices[x][1]);
      fade -= 1.0 / 3.0;
      //glColor4f(0, 1, 0,fade);
      glColor4f(color[0], color[1], color[2], fade);
    }
    //glColor4f(0,0,0,1);
    glVertex2f(_vertices[0][0], _vertices[0][1]);
  glEnd();

  glColor4f(0,0,1,1);
  glBegin(GL_POINTS);
    for (int x = 0; x < 3; x++)
      glVertex2f(_vertices[x][0], _vertices[x][1]);
  glEnd();
}
*/

///////////////////////////////////////////////////////////////////////
// draw a line using triangles
///////////////////////////////////////////////////////////////////////
void drawLine(const VEC2& v0, const VEC2& v1, const float lineWidth = 0.0025)
{
  glLineWidth(1);
  glBegin(GL_LINES);
    glVertex3f(v0[0], v0[1], 1.0);
    glVertex3f(v1[0], v1[1], 1.0);
  glEnd();
  /*
  static GLUquadricObj* quadric = gluNewQuadric();

  // endcaps
  glPushMatrix();
    glTranslatef(v0[0], v0[1], 0.0); 
    gluDisk(quadric, 0.0, lineWidth * 0.5, 20, 1);
  glPopMatrix();
  
  glPushMatrix();
    glTranslatef(v1[0], v1[1], 0.0); 
    gluDisk(quadric, 0.0, lineWidth * 0.5, 20, 1);
  glPopMatrix();

  VEC2 diff = v1 - v0;
  VEC2 translate = v0;
  if (v0[0] > v1[0])
  {
    diff = v0 - v1;
    translate = v1;
  }

  const float angle = (atan(diff[1] / diff[0]) / (2.0 * M_PI)) * 360.0;
  const float length = diff.norm();

  // bumped up vertices
  const VEC2 zero(0,    -lineWidth * 0.5);
  const VEC2 one(length,-lineWidth * 0.5);
  const VEC2 two(length, lineWidth * 0.5);
  const VEC2 three(0,    lineWidth * 0.5);

  glPushMatrix();
    glTranslatef(translate[0], translate[1], 0.1);
    glRotatef(angle, 0,0,1.0); 
    glBegin(GL_TRIANGLES);
      glVertex2f(zero[0], zero[1]);
      glVertex2f(one[0], one[1]);
      glVertex2f(two[0], two[1]);
    
      glVertex2f(two[0], two[1]);
      glVertex2f(three[0], three[1]);
      glVertex2f(zero[0], zero[1]);
    glEnd();
  glPopMatrix();
  */
}

///////////////////////////////////////////////////////////////////////
// Draw the square to GL
///////////////////////////////////////////////////////////////////////
void TRIANGLE::draw(const VEC4& color)
{
  MATRIX2 F = computeF();
  Real J = F.determinant();
  bool inverted = (J < 0.0) ? true : false;

  if (!inverted)
    glColor4f(color[0], color[1], color[2], color[3]);
  else
    glColor4f(254.0 / 255, 140.0 / 255, 117.0 / 255, 1);
  glBegin(GL_TRIANGLES);
    glVertex2f((*_vertices[0])[0], (*_vertices[0])[1]);
    glVertex2f((*_vertices[1])[0], (*_vertices[1])[1]);
    glVertex2f((*_vertices[2])[0], (*_vertices[2])[1]);
  glEnd();
  
  glColor4f(0,0,0,1);
  const Real lineWidth = 0.0025;
  drawLine(*_vertices[0], *_vertices[1], lineWidth);
  drawLine(*_vertices[1], *_vertices[2], lineWidth);
  drawLine(*_vertices[2], *_vertices[0], lineWidth);
}

///////////////////////////////////////////////////////////////////////
// Draw the square to GL
///////////////////////////////////////////////////////////////////////
void TRIANGLE::drawOutline(const Real lineWidth)
{
  // try drawing with quadrics and triangles instead
  glColor4f(0,0,0,1);
  //const Real lineWidth = 0.0025;
  drawLine(*_vertices[0], *_vertices[1], lineWidth);
  drawLine(*_vertices[1], *_vertices[2], lineWidth);
  drawLine(*_vertices[2], *_vertices[0], lineWidth);
}

/*
///////////////////////////////////////////////////////////////////////
// Draw the square forces to GL
///////////////////////////////////////////////////////////////////////
void TRIANGLE::drawForces()
{
  glLineWidth(2);
  glBegin(GL_LINES);
    for (int x = 0; x < 3; x++)
    {
      glColor4f(1,0,0,1);
      glVertex2f(_vertices[x][0], _vertices[x][1]);

      Real dt = 0.001;
      VEC2 forceEnd = _vertices[x] + _forces.col(x) * dt;
      
      glColor4f(1,0,0,0);
      glVertex2f(forceEnd[0], forceEnd[1]);
    }
  glEnd();
}
*/

/*
///////////////////////////////////////////////////////////////////////
// get the shape function derivative matrix H for computing F
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::computeH(const VEC2& Xi)
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
*/

///////////////////////////////////////////////////////////////////////
// build the material coordinate matrix
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::vertexMatrix() const
{
  /*
  MATRIX final(2,2);
  final(0,0) = _vertices[0][0];
  final(0,1) = _vertices[1][0];
  final(0,2) = _vertices[2][0];
  //final(0,3) = _vertices[3][0];
  
  final(1,0) = _vertices[0][1];
  final(1,1) = _vertices[1][1];
  final(1,2) = _vertices[2][1];
  //final(1,3) = _vertices[3][1];
  */
  MATRIX final(2,2);
  final.col(0) = (*_vertices[1]) - (*_vertices[0]);
  final.col(1) = (*_vertices[2]) - (*_vertices[0]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the deformation gradient at Xi, the normalized position
// inside a [-1,-1], [1,1] square
///////////////////////////////////////////////////////////////////////
MATRIX2 TRIANGLE::computeF() const
{
  MATRIX Ds = vertexMatrix();
  /*
  MATRIX H = computeH(Xi);

  MATRIX2 DmH_Inv = (_Dm * H).inverse();
  */

  //return (Ds * H) * DmH_Inv;
  return Ds * _Dm.inverse();
}

///////////////////////////////////////////////////////////////////////
// compute the strain energy
///////////////////////////////////////////////////////////////////////
Real TRIANGLE::psi() const
{
  MATRIX2 F = computeF();

  return _material->psi(F);
}

///////////////////////////////////////////////////////////////////////
// scale whole square
///////////////////////////////////////////////////////////////////////
void TRIANGLE::scale(const VEC2& scalar)
{
  for (int x = 0; x < 3; x++)
  {
    _vertices[x][0] *= scalar[0];
    _vertices[x][1] *= scalar[1];
  }
}

///////////////////////////////////////////////////////////////////////
// translate the whole square
///////////////////////////////////////////////////////////////////////
void TRIANGLE::translate(const VEC2& translation)
{
  for (int x = 0; x < 3; x++)
    (*_vertices[x]) += translation;
}

///////////////////////////////////////////////////////////////////////
// rotate the whole square
///////////////////////////////////////////////////////////////////////
void TRIANGLE::rotate(const Real theta)
{
  VEC2 center = vertexAverage();

  MATRIX2 rotation;
  rotation(0,0) = rotation(1,1) = cos(theta);
  rotation(0,1) = -sin(theta);
  rotation(1,0) = sin(theta);

  for (int x = 0; x < 3; x++)
    (*_vertices[x]) = rotation * ((*_vertices[x]) - center) + center;
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

/*
///////////////////////////////////////////////////////////////////////
// get the average F at all the Gauss points
///////////////////////////////////////////////////////////////////////
MATRIX2 TRIANGLE::gaussF()
{
  const Real invSqrt3 = 1.0 / sqrt(3.0);
  const VEC2 gaussPoints[] = {invSqrt3 * _vertices[0],
                              invSqrt3 * _vertices[1],
                              invSqrt3 * _vertices[2],
                              invSqrt3 * _vertices[3]};

  MATRIX2 final = computeF(gaussPoints[0]);
  for (int x = 1; x < 4; x++)
    final += computeF(gaussPoints[x]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// populate the force vector
///////////////////////////////////////////////////////////////////////
void TRIANGLE::computeForces()
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
*/

///////////////////////////////////////////////////////////////////////
// Compure Moore-Penrose pseudo-inverse
///////////////////////////////////////////////////////////////////////
static MATRIX2 pinv(const MATRIX2& A)
{
  // get the degenerate index
  JacobiSVD<MATRIX2> svd(A, ComputeFullU | ComputeFullV);
  VECTOR sigmas = svd.singularValues();
  //cout << " Singulars: " << sigmas.transpose() << endl;
  for (int x = 0; x < sigmas.size(); x++)
  {
    if (fabs(sigmas[x]) > 1e-8)
      sigmas[x] = 1.0 / sigmas[x];
    else
      sigmas[x] = 0.0;
  }

  return svd.matrixV() * sigmas.asDiagonal() * svd.matrixU().transpose();
}

///////////////////////////////////////////////////////////////////////
// find out which edge is least degenerate edge
///////////////////////////////////////////////////////////////////////
int TRIANGLE::findLeastDegenerateEdge() const
{
  // get the lengths
  Real lengths[3];
  for (int x = 0; x < 3; x++)
  {
    VEC2 diff = _restPose[(x + 1) % 3] - _restPose[x];
    lengths[x] = diff.norm();
  }

  // find the longest one
  int bestEdge = 0;
  Real bestLength = lengths[0];
  for (int x = 1; x < 3; x++)
    if (lengths[x] > bestLength)  // this will bias towards picking edge 1
    //if (lengths[x] >= bestLength) // this will bias towards picking edge 2
    {
      bestEdge = x;
      bestLength = lengths[x];
    }

#if 0
  cout << " Best found edge: " << bestEdge << endl;

  if (bestEdge == 1)
  {
    for (int x = 0; x < 3; x++)
      cout << " rest " << x << ": " << _restPose[x] << endl;
  }
#endif

  return bestEdge;
}

///////////////////////////////////////////////////////////////////////
// peek at the Dm and some alternatives
///////////////////////////////////////////////////////////////////////
void TRIANGLE::analyzeDm(const MATRIX2& Dm) const
{
  JacobiSVD<MATRIX2> svdDm(Dm, ComputeFullU | ComputeFullV);
  cout << endl;
  cout << "========================================== " << endl;
  cout << " Original Dm: " << endl;
  cout << Dm << endl;
  /*
  cout << " Original SVD: " << endl;
  cout << " U: " << endl << svdDm.matrixU() << endl;
  cout << " V: " << endl << svdDm.matrixV() << endl;
  cout << " sigmas: " << endl << svdDm.singularValues() << endl;
  cout << " Dm: " << endl << Dm << endl;
  */

  const int goodEdge = findLeastDegenerateEdge();
  int zero = goodEdge;
  int one = (zero + 1) % 3;
  int two = (zero + 2) % 3;   
  VEC2 edgeBar0 = _restPose[one] - _restPose[zero];
  VEC2 edgeBar2 = _restPose[two] - _restPose[zero];
  MATRIX2 DmOriginal;
  DmOriginal.col(0) = edgeBar0;
  DmOriginal.col(1) = edgeBar2;
  /*
  Real coordinate = edgeBar2.dot(edgeBar0) / edgeBar0.norm();
  VEC2 mBar = coordinate * edgeBar0.normalized();
  VEC2 cBar = _restPose[one] - _restPose[zero];

  // this just does a M_PI / 2 rotation of cBar
  VEC2 dBar = VEC2(-cBar[1], cBar[0]);
  VEC2 dBarNormalized = dBar / dBar.norm();
  //double eTall = sqrt(3.0) * 0.5 * cBar.norm();
  double eTall = cBar.norm();
  VEC2 phantomBar = mBar + eTall * dBarNormalized;

  VEC2 DmCol = phantomBar;
  MATRIX2 DmIsotropic;
  DmIsotropic.col(0) = cBar;
  DmIsotropic.col(1) = DmCol;
  cout << " Current Dm: " << endl;
  cout << DmIsotropic << endl;
  */

  MATRIX2 DmCopy = Dm; 
  DmCopy.col(1).setZero();

  VEC2 cBar = Dm.col(0);
  //VEC2 cn = cBar.normalized();
  Real cMag = cBar.norm();

  MATRIX2 T;
  T << 0, -1, 1,  0;

  //VEC2 col2 = (1.0 / (cMag * cMag)) * (cn * cn.transpose()) * edgeBar2 + T * cBar;
  VEC2 col2 = (1.0 / (cMag * cMag)) * (cBar * cBar.transpose()) * edgeBar2 + T * cBar;
  DmCopy.col(1) = col2;
  cout << " Dm copy: " << endl << DmCopy << endl;
  cout << endl;
  cout << " c magnitude: " << cMag << endl;
  VEC2 after = (1.0 / (cMag * cMag)) * (cBar * cBar.transpose()) * cBar;
  cout << " before: " << cBar.transpose() << endl;
  cout << " after:  " << after.transpose() << endl;

  MATRIX2 C = (1.0 / (cMag * cMag)) * (cBar * cBar.transpose());

  //MATRIX2 I = MATRIX2::Identity();
  MATRIX2 E;
  E << 0, 1.0, 0,0;
  MATRIX2 rightMatrix = T * DmOriginal * E;
  MATRIX2 DmAgain = C * DmOriginal + rightMatrix;
  cout << " Outer, after product: " << endl << C * DmOriginal << endl;
  cout << " Dm, again: " << endl << DmAgain << endl;

  /*
  JacobiSVD<MATRIX2> svdDmCopy(DmCopy, ComputeFullU | ComputeFullV);
  cout << endl << " Test SVD: " << endl;
  cout << " U: " << endl << svdDmCopy.matrixU() << endl;
  cout << " V: " << endl << svdDmCopy.matrixV() << endl;
  cout << " sigmas: " << endl << svdDmCopy.singularValues() << endl;
  cout << " Dm copy: " << endl << DmCopy << endl;
  */
  MATRIX2 diff = DmCopy - DmAgain;
  Real diffNorm = diff.norm();
  if (diffNorm > 1e-8)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " New computation failed! Diff magnitude: " << diffNorm << endl;
    cout << diff << endl;
    exit(0);
  }
  cout << "========================================== " << endl;


}

///////////////////////////////////////////////////////////////////////
// peek at the Ds and some alternatives
///////////////////////////////////////////////////////////////////////
void TRIANGLE::analyzeDs(const MATRIX2& Ds) const
{
  JacobiSVD<MATRIX2> svdDs(Ds, ComputeFullU | ComputeFullV);
  cout << endl;
  cout << "========================================== " << endl;
  cout << " Original Ds: " << endl;
  cout << Ds << endl;
  
  const int goodEdge = findLeastDegenerateEdge();
  int zero = goodEdge;
  int one = (zero + 1) % 3;
  int two = (zero + 2) % 3;   
  VEC2 edge0 = *_vertices[one] - *_vertices[zero];
  VEC2 edge2 = *_vertices[two] - *_vertices[zero];

  VEC2 c = edge0;
  VEC2 cBar = _restPose[one] - _restPose[zero];
  VEC2 eBar = _restPose[two] - _restPose[zero];
  MATRIX2 DmOriginal;
  DmOriginal.col(0) = cBar;
  DmOriginal.col(1) = eBar;
  MATRIX2 DsOriginal;
  DsOriginal.col(0) = edge0;
  DsOriginal.col(1) = edge2;

  // build some verification matrices
  MATRIX2 T;
  T << 0, -1, 1,  0;
  MATRIX2 E0;
  E0 << 1.0, 0.0, 0,0;
  MATRIX2 E1;
  E1 << 0, 1,0,0;
  
  MATRIX2 Delta0;
  Delta0 << 1,0,0,0;
  MATRIX2 Delta1;
  Delta1 << 0,0,0,1;

  // recompute from scratch
  MATRIX2 DsRecomputed;
  VEC2 edgeBar2 = _restPose[two] - _restPose[zero];
  VEC2 phantom = ((edgeBar2.dot(c)) / c.squaredNorm()) * c + (cBar.norm() / c.norm()) * T * c;
  DsRecomputed.col(0).setZero();
  DsRecomputed.col(1) = phantom;
  DsRecomputed += DsOriginal * E0;

  MATRIX2 C = (c * c.transpose()) / c.squaredNorm();
  MATRIX2 DsHat = DsOriginal;
  DsHat.col(1) = eBar;
  //MATRIX2 DsGuess = C * DsHat + (cBar.norm() / c.norm()) * T * DsOriginal * E1;
  MATRIX2 DsGuess = C * (DsOriginal * Delta0 + DmOriginal * Delta1) + (cBar.norm() / c.norm())  * T * DsOriginal * E1;
  cout << " Ds guess: " << endl;
  cout << DsGuess << endl;

  MATRIX2 diff = Ds - DsGuess;
  Real diffNorm = diff.norm();
  if (diffNorm > 1e-8)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " New computation failed! Diff magnitude: " << diffNorm << endl;
    cout << diff << endl;

    cout << " Diff wrt recomputed: " << endl;
    cout << (Ds - DsRecomputed).norm() << endl;
    cout << " Diff wrt recomputed norm: " << endl;
    cout << (Ds - DsRecomputed) << endl;

    exit(0);
  }
  cout << "========================================== " << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE::computeDegenerateQuantities(MATRIX2& DmIsotropic,
                                           MATRIX2& DsIsotropic, 
                                           MATRIX2& DsAnisotropic, 
                                           VEC2& cBar, 
                                           Real& phantomRestArea,
                                           MATERIAL* materialPointer,
                                           MATRIX6& permutation, 
                                           bool debug) const
{
#if USING_ANISOTROPIC_ARAP
  ANISOTROPIC_ARAP& material = *((ANISOTROPIC_ARAP*)materialPointer);
#else
  ANISOTROPIC_DIRICHLET& material = *((ANISOTROPIC_DIRICHLET*)materialPointer);
#endif

  // depending on which edge is the "good" edge, recompute the indices
  const int goodEdge = findLeastDegenerateEdge();
  int zero = goodEdge;
  int one = (zero + 1) % 3;
  int two = (zero + 2) % 3;   // this vertex corresponds to the one in
                              // the degenerate direction
 
#define USING_NEW_PROJECTION 1

#if USING_NEW_PROJECTION
  VEC2 edgeBar0 = _restPose[one] - _restPose[zero];
  VEC2 edgeBar2 = _restPose[two] - _restPose[zero];
  Real coordinate = edgeBar2.dot(edgeBar0) / edgeBar0.norm();
  VEC2 mBar = _restPose[zero] + coordinate * edgeBar0.normalized();
#else
  VEC2 mBar = 0.5 * (_restPose[one] + _restPose[zero]);
#endif

  cBar      = _restPose[one] - _restPose[zero];

  // this just does a M_PI / 2 rotation of cBar
  VEC2 dBar = VEC2(-cBar[1], cBar[0]);
  VEC2 dBarNormalized = dBar / dBar.norm();
  //double eTall = sqrt(3.0) * 0.5 * cBar.norm();
  double eTall = cBar.norm();
  VEC2 phantomBar = mBar + eTall * dBarNormalized;

  VEC2 DmCol = phantomBar - _restPose[zero];
  DmIsotropic.col(0) = cBar;

// show what happens when you use just an orthogonal frame
#define USING_D_COLUMN 0

#if USING_D_COLUMN
  DmIsotropic.col(1) = dBarNormalized;
#else
  DmIsotropic.col(1) = DmCol;
#endif

  static bool firstTime = true;

  if (firstTime)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " DM CHECK DEACTIVATED " << endl;
    firstTime = false;
  }
  //analyzeDm(DmIsotropic);

#if USING_ANISOTROPIC_ARAP
  VEC3 restPose3[3];
  for (int x = 0; x < 3; x++)
  {
    restPose3[x].setZero();
    for (int y = 0; y < 2; y++)
      restPose3[x][y] = _restPose[x][y];
  }
  VEC3 triangleNormal = (restPose3[2] - restPose3[0]).cross(restPose3[1] - restPose3[0]);
  phantomRestArea = triangleNormal.norm() * 0.5;
#else
  // getting the phantom rest area
  // first, get triangle normal in R^3
  VEC3 triangleNormal;
  VEC3 restPose3[3];
  for (int x = 0; x < 3; x++)
  {
    restPose3[x].setZero();
    for (int y = 0; y < 2; y++)
      restPose3[x][y] = _restPose[x][y];
  }

  restPose3[two][0] = phantomBar[0];
  restPose3[two][1] = phantomBar[1];

  // this does not need to be made generic, because
  // the area formula remains the same
  triangleNormal = (restPose3[2] - restPose3[0]).cross(restPose3[1] - restPose3[0]);
  phantomRestArea = triangleNormal.norm() * 0.5;
#endif
  
  VEC2 c =        *_vertices[one] - *_vertices[zero];
#if USING_NEW_PROJECTION
  VEC2 edge0 = *_vertices[one] - *_vertices[zero];

  // normalization doesn't seem to matter here -- all that does
  // is the position of the phantom rest vertex. It'll always try to 
  // pull back there.
  VEC2 m = *_vertices[zero] + coordinate * edge0.normalized();
#else
  VEC2 m = 0.5 * (*_vertices[one] + *_vertices[zero]);
#endif

  // this just does a M_PI / 2 rotation of c
  VEC2 d = VEC2(-c[1], c[0]);
  VEC2 dNormalized = d / d.norm();
  VEC2 phantom = (m + eTall * dNormalized);

  VEC2 DsCol = phantom - *_vertices[zero];
  DsIsotropic.col(0) = c;
#if USING_D_COLUMN
  DsIsotropic.col(1) = dNormalized;
#else
  DsIsotropic.col(1) = DsCol;
#endif

  if (firstTime)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " DS CHECK DEACTIVATED " << endl;
  }
  //analyzeDs(DsIsotropic);

  DsAnisotropic.col(0) = *_vertices[one] - *_vertices[zero];
  DsAnisotropic.col(1) = *_vertices[two] - *_vertices[zero];

  material.setAnisotropyDirection(dBarNormalized);

  //if (debug)
  if (0)
  {
    cout << " cBar: " << endl << cBar << endl;
    cout << " dBar: " << endl << dBar << endl;
    cout << " DmIsotropic: " << endl << DmIsotropic << endl;
    cout << " phantom rest:" << endl << phantomBar << endl;
#if USING_NEW_PROJECTION
    cout << " coordinate: " << endl << coordinate << endl;
#endif
    cout << " mBar:       " << endl << mBar << endl;
    cout << " mBar guess: " << endl << 0.5 * (_restPose[one] + _restPose[zero]) << endl;
    cout << " Phantom rest area: " << endl << phantomRestArea << endl;
    
    cout << " phantom: " << endl << phantom << endl;
    cout << " DsIsotropic: " << endl << DsIsotropic << endl;
    cout << " DsAnisotropic: " << endl << DsAnisotropic << endl;
    cout << " fiber direction: " << endl << dBarNormalized << endl;
  }

  // initialize assuming the good edge is 0
  permutation = MATRIX6::Identity();
  if (goodEdge == 1)
  {
    permutation.setZero();
    permutation(2,0) = 1.0;
    permutation(3,1) = 1.0;

    permutation(4,2) = 1.0;
    permutation(5,3) = 1.0;
    
    permutation(0,4) = 1.0;
    permutation(1,5) = 1.0;
  }
  if (goodEdge == 2)
  {
    permutation.setZero();
    permutation(4,0) = 1.0;
    permutation(5,1) = 1.0;

    permutation(0,2) = 1.0;
    permutation(1,3) = 1.0;
    
    permutation(2,4) = 1.0;
    permutation(3,5) = 1.0;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::computeAnisotropicDegenerateForce()
{
  MATRIX2 DmIsotropic;
  MATRIX2 DsIsotropic;
  MATRIX2 DsAnisotropic;
  VEC2 cBar;
  Real phantomRestArea;
#if USING_ANISOTROPIC_ARAP
  //ANISOTROPIC_ARAP material(SUPER_GLUE, SUPER_GLUE);
  ANISOTROPIC_ARAP material(superGlue, superGlue);
#else
  //ANISOTROPIC_DIRICHLET material(SUPER_GLUE, SUPER_GLUE);
  ANISOTROPIC_DIRICHLET material(superGlue, superGlue);
#endif
  MATRIX6 permutation;
  computeDegenerateQuantities(DmIsotropic, 
                              DsIsotropic,
                              DsAnisotropic,
                              cBar, 
                              phantomRestArea, 
                              &material,
                              permutation,
                              true);
  
  // anisotropic force stuff starts here
#if USING_ANISOTROPIC_ARAP
  const int goodEdge = findLeastDegenerateEdge();
  int zero = goodEdge;
  int one = (zero + 1) % 3;
  int two = (zero + 2) % 3;
  MATRIX2 DmAnisotropic;
  DmAnisotropic.col(0) = _restPose[one] - _restPose[zero];
  DmAnisotropic.col(1) = _restPose[two] - _restPose[zero];
#else
  MATRIX2 DmAnisotropic = DmIsotropic;
#endif
  
  MATRIX2 FAnisotropic = DsAnisotropic * DmAnisotropic.inverse();
  MATRIX2 PAnisotropic = material.PK1(FAnisotropic);

  MATRIX2 DmInvAnisotropic = DmAnisotropic.inverse();
  /*
  MATRIX2 DmPreinverse = DmAnisotropic;
  //DmPreinverse.col(0).setZero();
  MATRIX2 DmInvAnisotropic = pinv(DmPreinverse);
  */
  MATRIX pFpuTranspose  = permutation * pFpuDegenerate(DmInvAnisotropic).transpose();
  VECTOR forceAnisotopic = pFpuTranspose * flatten(PAnisotropic);

  /*
  // trial
  {
    MATRIX2 DmPreinverseTrial = DmAnisotropic;
    DmPreinverseTrial.col(0).setZero();
    MATRIX2 DmInvAnisotropicTrial = pinv(DmPreinverseTrial);
    MATRIX pFpuTransposeTrial  = permutation * pFpuDegenerate(DmInvAnisotropicTrial).transpose();
    VECTOR forceAnisotopicTrial = pFpuTransposeTrial * flatten(PAnisotropic);

    cout << " force:       " << forceAnisotopic.transpose() << endl;
    cout << " force trial: " << forceAnisotopicTrial.transpose() << endl;
  }
  */

  /*
  int edge = findLeastDegenerateEdge();
  if (edge == 0 || edge == 1)
  {
    cout << " Edge: " << edge  << endl;
    cout << " Anisotropy direction: " << endl << material.anisotropyDirection() << endl;
    cout << " Ds anisotropic: " << endl << DsAnisotropic << endl;
    cout << " Dm anisotropic: " << endl << DmAnisotropic << endl;
    cout << " F Anisotropic: " << endl << FAnisotropic << endl;
    JacobiSVD<MATRIX2> svd(FAnisotropic, ComputeFullU | ComputeFullV);
    cout << " Stretch: " << endl << svd.singularValues() << endl;
    cout << " P Anisotropic: " << endl << PAnisotropic << endl;
    cout << " force anisotropic: " << endl << forceAnisotopic << endl;
    cout << " Phantom Rest Area: " << endl << phantomRestArea << endl;
  }
  */

  // DEBUG: why does removing the minus stabilize things?
  return -1.0 * phantomRestArea * forceAnisotopic;
  //return phantomRestArea * forceAnisotopic;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::computeIsotropicDegenerateForce()
{
  MATRIX2 DmIsotropic;
  MATRIX2 DsIsotropic;
  MATRIX2 DsAnisotropic;
  VEC2 cBar;
  Real phantomRestArea;
#if USING_ANISOTROPIC_ARAP
  //ANISOTROPIC_ARAP anisotropic(SUPER_GLUE, SUPER_GLUE);
  ANISOTROPIC_ARAP anisotropic(superGlue, superGlue);
#else
  //ANISOTROPIC_DIRICHLET anisotropic(SUPER_GLUE, SUPER_GLUE);
  ANISOTROPIC_DIRICHLET anisotropic(superGlue, superGlue);
#endif
  MATRIX6 permutation;
  computeDegenerateQuantities(DmIsotropic, 
                              DsIsotropic,
                              DsAnisotropic,
                              cBar, 
                              phantomRestArea, 
                              &anisotropic,
                              permutation);
  
  MATRIX2 FIsotropic = DsIsotropic * DmIsotropic.inverse();
  MATRIX2 PIsotropic = _material->PK1(FIsotropic);

#if USING_DS_ZEROING
  MATRIX2 DmInvIsotropic = DmIsotropic.inverse();
  //MATRIX pFpuTranspose  = permutation * pFpuDegenerateIsotropic(DmInvIsotropic).transpose();
  MATRIX pFpuTranspose  = permutation * pFpuDegenerateComplete(DmInvIsotropic).transpose();

  //cout << " pFpuDegenerateComplete: " << endl << pFpuDegenerateComplete(DmInvIsotropic) << endl;
  //cout << " pFpuDegenerateIsotropic: " << endl << pFpuDegenerateIsotropic(DmInvIsotropic) << endl;
  //exit(0);
#else
  MATRIX2 DmPreinverse;
  DmPreinverse.setZero();
  DmPreinverse.col(0) = cBar;
  MATRIX2 DmInvIsotropic = pinv(DmPreinverse);
  
  MATRIX pFpuTranspose  = permutation * pFpuDegenerate(DmInvIsotropic).transpose();
#endif
  VECTOR forceIsotropic = pFpuTranspose * flatten(PIsotropic);
  return -1.0 * restArea() * forceIsotropic;
}

///////////////////////////////////////////////////////////////////////
// populate the force vector
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::computeForceVectorForDegenerateElement()
{
  return computeIsotropicDegenerateForce() + computeAnisotropicDegenerateForce();
  //return computeIsotropicDegenerateForce();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real TRIANGLE::computeAnisotropicDegeneratePsi() const
{
  MATRIX2 DmIsotropic;
  MATRIX2 DsIsotropic;
  MATRIX2 DsAnisotropic;
  VEC2 cBar;
  Real phantomRestArea;
#if USING_ANISOTROPIC_ARAP
  //ANISOTROPIC_ARAP material(SUPER_GLUE, SUPER_GLUE);
  ANISOTROPIC_ARAP material(superGlue, superGlue);
#else
  //ANISOTROPIC_DIRICHLET material(SUPER_GLUE, SUPER_GLUE);
  ANISOTROPIC_DIRICHLET material(superGlue, superGlue);
#endif
  MATRIX6 permutation;
  computeDegenerateQuantities(DmIsotropic, 
                              DsIsotropic,
                              DsAnisotropic,
                              cBar, 
                              phantomRestArea, 
                              &material,
                              permutation);

  // anisotropic force stuff starts here
#if USING_ANISOTROPIC_ARAP
  const int goodEdge = findLeastDegenerateEdge();
  int zero = goodEdge;
  int one = (zero + 1) % 3;
  int two = (zero + 2) % 3;
  MATRIX2 DmAnisotropic;
  DmAnisotropic.col(0) = _restPose[one] - _restPose[zero];
  DmAnisotropic.col(1) = _restPose[two] - _restPose[zero];
#else
  MATRIX2 DmAnisotropic = DmIsotropic;
#endif

  MATRIX2 FAnisotropic = DsAnisotropic * DmAnisotropic.inverse();
  Real psiAnisotropic = material.psi(FAnisotropic);
 
  return phantomRestArea * psiAnisotropic;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real TRIANGLE::computeIsotropicDegeneratePsi() const
{
  MATRIX2 DmIsotropic;
  MATRIX2 DsIsotropic;
  MATRIX2 DsAnisotropic;
  VEC2 cBar;
  Real phantomRestArea;
#if USING_ANISOTROPIC_ARAP
  //ANISOTROPIC_ARAP material(SUPER_GLUE, SUPER_GLUE);
  ANISOTROPIC_ARAP material(superGlue, superGlue);
#else
  //ANISOTROPIC_DIRICHLET material(SUPER_GLUE, SUPER_GLUE);
  ANISOTROPIC_DIRICHLET material(superGlue, superGlue);
#endif
  MATRIX6 permutation;
  computeDegenerateQuantities(DmIsotropic, 
                              DsIsotropic,
                              DsAnisotropic,
                              cBar, 
                              phantomRestArea, 
                              &material,
                              permutation);
  
  MATRIX2 FIsotropic = DsIsotropic * DmIsotropic.inverse();
  Real psiIsotropic = _material->psi(FIsotropic);

  return restArea() * psiIsotropic;
}

///////////////////////////////////////////////////////////////////////
// populate the force vector
///////////////////////////////////////////////////////////////////////
Real TRIANGLE::psiDegenerate() const
{
  return computeIsotropicDegeneratePsi() + computeAnisotropicDegeneratePsi();
}

///////////////////////////////////////////////////////////////////////
// populate the force vector
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::computeForceVectorIben()
{
  VECTOR forceVector(6);
  forceVector.setZero();

  MATRIX2 F = computeF();
  MATRIX2 P = _material->PK1(F);
  VECTOR flat = flatten(P);

  forceVector = pFpuDegenerate().transpose() * flat;

  return -1.0 * forceVector * restArea();
}

///////////////////////////////////////////////////////////////////////
// populate the force vector
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::computeForceVector()
{
  VECTOR forceVector(6);
  forceVector.setZero();

  MATRIX2 F = computeF();
  MATRIX2 P = _material->PK1(F);
  VECTOR flat = flatten(P);

  forceVector = pFpu().transpose() * flat;

  return -1.0 * forceVector * restArea();
}

///////////////////////////////////////////////////////////////////////
// populate the force vector
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::computeForceVectorFast()
{
  VECTOR forceVector(6);
  forceVector.setZero();

  MATRIX2 F = computeF();
  MATRIX2 P = _material->PK1(F);

  for (int x = 0; x < 6; x++)
  {
    MATRIX2 partialFpartialu = pFpu(x);
    forceVector[x] = P.cwiseProduct(partialFpartialu).sum();
  }

  return -1.0 * forceVector * restArea();
}

///////////////////////////////////////////////////////////////////////
// populate the force vector
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::computeForceVectorDEC()
{
  MATRIX2 F = computeF();
  MATRIX2 P = _material->PK1(F);

  // get the normals
  VEC2 rotatedEdges[3];
  for (int x = 0; x < 3; x++)
  {
    int next = (x + 1) % 3;
    VEC2 diff = _restPose[x] - _restPose[next];

    rotatedEdges[x][0] = diff[1];
    rotatedEdges[x][1] = -diff[0];
  }

  // compute the final force
  MATRIX dFdu(6,4);
  dFdu.setZero();
  for (int x = 0; x < 3; x++)
  {
    int next = (x + 2) % 3;
    VEC2 e = (rotatedEdges[x] + rotatedEdges[next]) * 0.5;
    int row0 = 2 * x;
    int row1 = row0 + 1;
    dFdu(row0, 0) = e[0];
    dFdu(row1, 1) = e[0];
    dFdu(row0, 2) = e[1];
    dFdu(row1, 3) = e[1];
  }
  VECTOR forceVector = dFdu * flatten(P);

  return forceVector;
}

///////////////////////////////////////////////////////////////////////
// populate the force vector
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::computeForceVectorFVM()
{
  MATRIX2 F = computeF();
  MATRIX2 P = _material->PK1(F);

  // get the edge lengths
  Real lengths[3];
  for (int x = 0; x < 3; x++)
  {
    const int next = (x + 1) % 3;
    const VEC2 diff = _restPose[x] - _restPose[next];
    lengths[x] = diff.norm();
    //cout << " x: " << x << " next: " << next << endl;
    //cout << " diff: " << diff.transpose() << endl;
    //cout << " rest " << x << ":" << _restPose[x].transpose() << endl;
    //cout << " rest " << next << ":" << _restPose[next].transpose() << endl;
    //cout << " length " << x << ": " << lengths[x] << endl;
  }

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
  triangleNormal *= 1.0 / triangleNormal.norm();

  // get the normals
  VEC2 normals[3];
  for (int x = 0; x < 3; x++)
  {
    int next = (x + 1) % 3;
    const VEC2 diff = _restPose[x] - _restPose[next];
    VEC3 diff3;
    diff3[0] = diff[0];
    diff3[1] = diff[1];
    diff3[2] = 0;

    VEC3 normal = triangleNormal.cross(diff3);
    normal *= 1.0 / normal.norm();

    normals[x][0] = normal[0];
    normals[x][1] = normal[1];
    //cout << " normal " << x << ": " << normals[x] << endl;
  }

  // compute the final force
  VECTOR forceVector(6);
  forceVector.setZero();
  for (int x = 0; x < 3; x++)
  {
    int next = (x + 2) % 3;
    VEC2 force = P * (lengths[x] * normals[x] + lengths[next] * normals[next]);
    forceVector[2 * x] = force[0];
    forceVector[2 * x + 1] = force[1];
    //cout << " force " << x << ": " << force << endl;
  }

  //forceVector *= -1.0 / 3.0;
  forceVector *= 1.0 / 2.0;

  return forceVector;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::computeAnisotropicDegenerateForceJacobian()
{
  MATRIX2 DmIsotropic;
  MATRIX2 DsIsotropic;
  MATRIX2 DsAnisotropic;
  VEC2 cBar;
  Real phantomRestArea;
#if USING_ANISOTROPIC_ARAP
  //ANISOTROPIC_ARAP material(SUPER_GLUE, SUPER_GLUE);
  ANISOTROPIC_ARAP material(superGlue, superGlue);
#else
  //ANISOTROPIC_DIRICHLET material(SUPER_GLUE, SUPER_GLUE);
  ANISOTROPIC_DIRICHLET material(superGlue, superGlue);
#endif
  MATRIX6 permutation;
  computeDegenerateQuantities(DmIsotropic, 
                              DsIsotropic,
                              DsAnisotropic,
                              cBar, 
                              phantomRestArea, 
                              &material,
                              permutation);

  // anisotropic force stuff starts here
#if USING_ANISOTROPIC_ARAP
  const int goodEdge = findLeastDegenerateEdge();
  int zero = goodEdge;
  int one = (zero + 1) % 3;
  int two = (zero + 2) % 3;
  MATRIX2 DmAnisotropic;
  DmAnisotropic.col(0) = _restPose[one] - _restPose[zero];
  DmAnisotropic.col(1) = _restPose[two] - _restPose[zero];
#else
  MATRIX2 DmAnisotropic = DmIsotropic;
#endif

  MATRIX2 FAnisotropic = DsAnisotropic * DmAnisotropic.inverse();
  MATRIX2 DmInvAnisotropic = DmAnisotropic.inverse();
  MATRIX dpdf = material.DPDF(FAnisotropic);
  /*
  int edge = findLeastDegenerateEdge();
  if (edge == 1)
  {
    cout << " DPDF: " << endl << dpdf << endl;
  }
  */
 
  MATRIX pFpuTranspose = permutation * pFpuDegenerate(DmInvAnisotropic).transpose();
  MATRIX hessian = pFpuTranspose * dpdf * pFpuTranspose.transpose();

  return -1.0 * phantomRestArea * hessian;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::computeIsotropicDegenerateForceJacobian()
{
  MATRIX2 DmIsotropic;
  MATRIX2 DsIsotropic;
  MATRIX2 DsAnisotropic;
  VEC2 cBar;
  Real phantomRestArea;
#if USING_ANISOTROPIC_ARAP
  //ANISOTROPIC_ARAP material(SUPER_GLUE, SUPER_GLUE);
  ANISOTROPIC_ARAP material(superGlue, superGlue);
#else
  //ANISOTROPIC_DIRICHLET material(SUPER_GLUE, SUPER_GLUE);
  ANISOTROPIC_DIRICHLET material(superGlue, superGlue);
#endif
  MATRIX6 permutation;
  computeDegenerateQuantities(DmIsotropic, 
                              DsIsotropic,
                              DsAnisotropic,
                              cBar, 
                              phantomRestArea, 
                              &material,
                              permutation);

  MATRIX2 FIsotropic = DsIsotropic * DmIsotropic.inverse();
  MATRIX dpdf = _material->DPDF(FIsotropic);

#if USING_DS_ZEROING
  MATRIX2 DmInvIsotropic = DmIsotropic.inverse();
  //MATRIX pFpuTranspose = permutation * pFpuDegenerateIsotropic(DmInvIsotropic).transpose();
  MATRIX pFpuTranspose = permutation * pFpuDegenerateComplete(DmInvIsotropic).transpose();
#else
  MATRIX2 DmPreinverse;
  DmPreinverse.setZero();
  DmPreinverse.col(0) = cBar;
  MATRIX2 DmInvIsotropic = pinv(DmPreinverse);
 
  MATRIX pFpuTranspose = permutation * pFpuDegenerate(DmInvIsotropic).transpose();
#endif

  MATRIX hessian = pFpuTranspose * dpdf * pFpuTranspose.transpose();
  return -1.0 * restArea() * hessian;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::computeForceJacobianForDegenerateElement()
{
  return computeIsotropicDegenerateForceJacobian() +
         computeAnisotropicDegenerateForceJacobian();
  //return computeIsotropicDegenerateForceJacobian();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::computeForceJacobianIben()
{
  MATRIX2 F = computeF();
  MATRIX DPDF = _material->DPDF(F);
  MATRIX DFDu = pFpuDegenerate();

  MATRIX final(6,6);
  for (int j = 0; j < 6; j++)
    for (int i = 0; i < 6; i++)
      final(i,j) = (DPDF * DFDu.col(i)).dot(DFDu.col(j));

  return -1.0 * restArea() * final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::computeForceJacobianFast()
{
  MATRIX2 F = computeF();
  MATRIX DPDF = _material->DPDF(F);
  MATRIX DFDu = pFpu();

  MATRIX final(6,6);
  for (int j = 0; j < 6; j++)
    for (int i = 0; i < 6; i++)
      final(i,j) = (DPDF * DFDu.col(i)).dot(DFDu.col(j));

  return -1.0 * restArea() * final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::computeForceJacobianDEC()
{
  MATRIX2 F = computeF();
  MATRIX DPDF = _material->DPDF(F);

  // get the normals
  VEC2 rotatedEdges[3];
  for (int x = 0; x < 3; x++)
  {
    int next = (x + 1) % 3;
    VEC2 diff = _restPose[x] - _restPose[next];

    rotatedEdges[x][0] = diff[1];
    rotatedEdges[x][1] = -diff[0];
  }

  // compute the final force
  MATRIX dFdu(6,4);
  dFdu.setZero();
  for (int x = 0; x < 3; x++)
  {
    int next = (x + 2) % 3;
    VEC2 e = (rotatedEdges[x] + rotatedEdges[next]) * 0.5;
    int row0 = 2 * x;
    int row1 = row0 + 1;
    dFdu(row0, 0) = e[0];
    dFdu(row1, 1) = e[0];
    dFdu(row0, 2) = e[1];
    dFdu(row1, 3) = e[1];
  }

  MATRIX hessian = dFdu * DPDF * dFdu.transpose();
  assert(hessian.rows() == 6);
  assert(hessian.cols() == 6);
  return -1.0 * hessian / restArea();
}


/*
///////////////////////////////////////////////////////////////////////
// get the forces at a Gauss point
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::computeForceVector(const VEC2& gaussPoint)
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
MATRIX TRIANGLE::computeForces(const VEC2& gaussPoint)
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
MATRIX TRIANGLE::computeDPDu()
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
    final.col(i) = TRIANGLE::flatten(col);
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the forces Jacobian at a Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::computeForceJacobian()
{
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
// get PK2 over all the Gauss points
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::secondPiolaKirchhoff()
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
MATRIX TRIANGLE::computeForceJacobian(const VEC2& gaussPoint)
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

  // do it per entry
  MATRIX K(8,8);
  for (int i = 0; i < 8; i++)
  {
    MATRIX DFDu_i = computeDFDu(gaussPoint, i);
    VECTOR flatDFDu_i = flatten(DFDu_i);

    VECTOR flatDPDu = A * flatDFDu_i;

    MATRIX DPDu = unflatten(flatDPDu, 2);
    MATRIX squareKi = DPDu * Bmg;

    K.col(i) = flatten(squareKi);
  }

  return -1.0 * K;
}

///////////////////////////////////////////////////////////////////////
// get the PK1 Jacobian at a Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::computeDPDu(const VEC2& gaussPoint, const int column)
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

  // do it per entry
  //MATRIX K(8,8);
  //for (int i = 0; i < 8; i++)
  MATRIX DFDu_i = computeDFDu(gaussPoint, column);
  VECTOR flatDFDu_i = flatten(DFDu_i);

  VECTOR flatDPDu = A * flatDFDu_i;

  MATRIX DPDu = unflatten(flatDPDu, 2);

  return DPDu;
}

///////////////////////////////////////////////////////////////////////
// get the PK2 Jacobian at a Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::computeDSDu(const VEC2& gaussPoint, const int column)
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
MATRIX TRIANGLE::firstPiolaKirchhoff(const VEC2& gaussPoint)
{
  MATRIX2 F = computeF(gaussPoint);
  return _material->PK1(F);
}

///////////////////////////////////////////////////////////////////////
// get the second Piola Kirchhoff for a single Gauss point
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::secondPiolaKirchhoff(const VEC2& gaussPoint)
{
  MATRIX2 F = computeF(gaussPoint);
  return _material->PK2(F);
}
*/

/*
///////////////////////////////////////////////////////////////////////
// get the first Piola Kirchhoff for St. Venant-Kirchhoff
//
// this doesn't just call secondPiolaKirchhoffStVK because it would
// then build F twice
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::firstPiolaKirchhoffStVK(const VEC2& gaussPoint)
{
  MATRIX2 F = computeF(gaussPoint);
  MATRIX2 E = 0.5 * (F.transpose() * F - MATRIX2::Identity());

  // settings from the 2008 cubature paper
  const Real lambda = 1000;
  const Real mu = 5000;

  MATRIX2 S = lambda * E.trace() * MATRIX2::Identity() + 2.0 * mu * E;
  return F * S;
}

///////////////////////////////////////////////////////////////////////
// get the second Piola Kirchhoff for St. Venant-Kirchhoff
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::secondPiolaKirchhoffStVK(const VEC2& gaussPoint)
{
  MATRIX2 F = computeF(gaussPoint);
  MATRIX2 E = 0.5 * (F.transpose() * F - MATRIX2::Identity());

  // settings from the 2008 cubature paper
  const Real lambda = 1000;
  const Real mu = 5000;

  MATRIX2 S = lambda * E.trace() * MATRIX2::Identity() + 2.0 * mu * E;
  return S;
}

///////////////////////////////////////////////////////////////////////
// get the first Piola Kirchhoff for Co-rotational
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::firstPiolaKirchhoffCorotational(const VEC2& gaussPoint)
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

/*
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::forceVector() const
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
VECTOR TRIANGLE::getDisplacement() const
{
  VECTOR final(8);

  int i = 0;
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 2; y++, i++)
      final[i] = _vertices[x][y] - _restPose[x][y];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE::setDisplacement(VECTOR& u)
{
  assert(u.size() == 8);
  int i = 0;
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 2; y++, i++)
      _vertices[x][y] = _restPose[x][y] + u[i];
}
*/

/*
///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::flatten(const MATRIX& A)
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
MATRIX TRIANGLE::unflatten(const VECTOR& v, const int rows)
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
MATRIX TRIANGLE::computeDFDu(const VEC2& Xi)
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
MATRIX TRIANGLE::computeDFDu(const VEC2& Xi, const int index)
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
// compute the derivative of E = 1/2 (F^T * F - I)
// with respect to F
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::computeDEDF(const MATRIX2& F)
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
*/

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE::applyToRest(const MATRIX2& transform)
{
  for (int x = 0; x < 3; x++)
    (*_vertices[x]) = transform * _restPose[x];
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::flatten(const MATRIX2& A)
{
  VECTOR final(4);
  final[0] = A(0,0);
  final[1] = A(1,0);
  final[2] = A(0,1);
  final[3] = A(1,1);
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::pFpuDegenerate()
{
  MATRIX final(4,6);

  MATRIX2 DmInverse = _Dm.inverse();
  JacobiSVD<MATRIX2> svd(DmInverse, ComputeFullU | ComputeFullV);
  VECTOR sigmas = svd.singularValues();
  int degenerateIndex = -1;
  for (int x = 0; x < sigmas.size(); x++)
  {
    if (sigmas[x] > 1e8)
      degenerateIndex = x;
  }
  sigmas[degenerateIndex] = 0;
  DmInverse = svd.matrixU() * sigmas.asDiagonal() * svd.matrixV().transpose();
 
  for (int x = 0; x < 6; x++)
  {
    MATRIX2 column;
    column.setZero();

    if (x == 0)
    {
      column(0,0) = -1;
      column(0,1) = -1;
    }
    if (x == 1)
    {
      column(1,0) = -1;
      column(1,1) = -1;
    }
    if (x == 2)
    {
      column(0,0) = 1;
    }
    if (x == 3)
    {
      column(1,0) = 1;
    }
    if (x == 4)
    {
      column(0,1) = 1;
    }
    if (x == 5)
      column(1,1) = 1;

    column = column * DmInverse;

    final.col(x) = flatten(column);
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::pFpuDegenerateIsotropic(const MATRIX2& DmInverse)
{
  MATRIX final(4,6);

  for (int x = 0; x < 6; x++)
  {
    MATRIX2 column;
    column.setZero();

    if (x == 0)
    {
      column(0,0) = -1;
      column(0,1) = -1;
    }
    if (x == 1)
    {
      column(1,0) = -1;
      column(1,1) = -1;
    }
    if (x == 2)
    {
      column(0,0) = 1;
    }
    if (x == 3)
    {
      column(1,0) = 1;
    }
    /*
    if (x == 4)
    {
      column(0,1) = 1;
    }
    if (x == 5)
    {
      column(1,1) = 1;
    }
    */

    column = column * DmInverse;

    final.col(x) = flatten(column);
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::pFpuDegenerate(const MATRIX2& DmInverse)
{
  MATRIX final(4,6);

  for (int x = 0; x < 6; x++)
  {
    MATRIX2 column;
    column.setZero();

    if (x == 0)
    {
      column(0,0) = -1;
      column(0,1) = -1;
    }
    if (x == 1)
    {
      column(1,0) = -1;
      column(1,1) = -1;
    }
    if (x == 2)
    {
      column(0,0) = 1;
    }
    if (x == 3)
    {
      column(1,0) = 1;
    }
    if (x == 4)
    {
      column(0,1) = 1;
    }
    if (x == 5)
    {
      column(1,1) = 1;
    }

    column = column * DmInverse;

    final.col(x) = flatten(column);
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::pFpu()
{
  MATRIX final(4,6);

  for (int x = 0; x < 6; x++)
  {
    MATRIX2 column = pFpu(x);
    final.col(x) = flatten(column);
  }

  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " DFDu: " << endl << final.transpose()  * final << endl;

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX2 TRIANGLE::pFpu(const int index)
{
  MATRIX2 final;
  final.setZero();

  if (index == 0)
  {
    final(0,0) = -1;
    final(0,1) = -1;
  }
  if (index == 1)
  {
    final(1,0) = -1;
    final(1,1) = -1;
  }
  if (index == 2)
  {
    final(0,0) = 1;
  }
  if (index == 3)
  {
    final(1,0) = 1;
  }
  if (index == 4)
  {
    final(0,1) = 1;
  }
  if (index == 5)
    final(1,1) = 1;

  final = final * _Dm.inverse();

  return final;
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
VECTOR TRIANGLE::computeForceVectorPK2()
{
  MATRIX2 F = computeF();
  MATRIX2 PK2 = _material->PK2(F);
  MATRIX3 alpha;
  alpha.setOnes();
  for (unsigned int x = 0; x < 3; x++)
    for (unsigned int y = 0; y < 2; y++)
      alpha(y,x) = _restPose[x][y];
  MATRIX3 alphaInv = alpha.inverse();
  //MATRIX3 beta = alpha.inverse();
  MATRIX beta(3,2);
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 2; y++)
      beta(x,y) = alphaInv(x,y);

  VECTOR p(6);  
  for (unsigned int i = 0; i < 3; i++)
  {
    VEC2 v = *(_vertices[i]);
    p[2 * i] = v[0];
    p[2 * i + 1] = v[1];
  }

  VECTOR forceVector(6);
  MATRIX meat(6,6);
  meat.setZero();
  for (unsigned int i = 0; i < 3; i++)
  {
    //VEC2 force;
    //force.setZero();

    meat.setZero();
    for (unsigned int j = 0; j < 3; j++)
    {
      Real coeff = 0.0;
      for (unsigned int k = 0; k < 2; k++)
        for (unsigned int l = 0; l < 2; l++)
          coeff += beta(j,l) * beta(i,k) * PK2(k,l);

      meat(2 * i, 2 * j) = coeff;
      meat(2 * i + 1, 2 * j + 1) = coeff;
    }
    //forceVector[2 * i] = force[0];
    //forceVector[2 * i + 1] = force[1];
  }
  forceVector = meat * p;
  forceVector *= restArea();
  forceVector *= -1.0;

  /*
  VECTOR forceVector(6);
  for (unsigned int i = 0; i < 3; i++)
  {
    VEC2 force;
    force.setZero();

    MATRIX meat(2,6);
    meat.setZero();
    for (unsigned int j = 0; j < 3; j++)
    {
      Real coeff = 0.0;
      for (unsigned int k = 0; k < 2; k++)
        for (unsigned int l = 0; l < 2; l++)
          coeff += beta(j,l) * beta(i,k) * PK2(k,l);

      meat(0, 2 * j) = coeff;
      meat(1, 2 * j + 1) = coeff;
    }
    cout << " meat: " << endl << meat << endl;
    force = meat * p;

    forceVector[2 * i] = force[0];
    forceVector[2 * i + 1] = force[1];
  }
  forceVector *= restArea();
  forceVector *= -1.0;
  */
  return forceVector;
}
/*
{
  MATRIX2 F = computeF();
  MATRIX2 PK2 = _material->PK2(F);
  MATRIX3 alpha;
  alpha.setOnes();
  for (unsigned int x = 0; x < 3; x++)
    for (unsigned int y = 0; y < 2; y++)
      alpha(y,x) = _restPose[x][y];
  MATRIX3 beta = alpha.inverse();

  MATRIX G(3,2);
  G.setZero();
  G(0,0) = 1.0;
  G(1,1) = 1.0;

  MATRIX meat = beta * G * PK2 * G.transpose() * beta.transpose();
  VECTOR forceVector(6);
  for (unsigned int i = 0; i < 3; i++)
  {
    //VEC2 p = *(_vertices[i]);
    VEC2 p = *(_vertices[i]) - _restPose[i];
    // homogeneous version
    VEC3 pH;
    for (int x = 0; x < 2; x++)
      pH[x] = p[x];
    pH[2] = 1.0;
    VEC3 forceH = meat * pH;

    forceVector[2 * i] = forceH[0];
    forceVector[2 * i + 1] = forceH[1];
  }
  forceVector *= restArea();
  forceVector *= -1.0;
  return forceVector;
}
*/
/*
{
  MATRIX2 F = computeF();
  MATRIX2 PK2 = _material->PK2(F);
  MATRIX3 alpha;
  alpha.setOnes();
  for (unsigned int x = 0; x < 3; x++)
    for (unsigned int y = 0; y < 2; y++)
      alpha(y,x) = _restPose[x][y];
  MATRIX3 alphaInv = alpha.inverse();
  //MATRIX3 beta = alpha.inverse();
  MATRIX beta(3,2);
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 2; y++)
      beta(x,y) = alphaInv(x,y);

  VECTOR forceVector(6);
  for (unsigned int i = 0; i < 3; i++)
  {
    VEC2 force;
    force.setZero();
    for (unsigned int j = 0; j < 3; j++)
    {
      VEC2 p = *(_vertices[j]);
      Real coeff = 0.0;
      for (unsigned int k = 0; k < 2; k++)
        for (unsigned int l = 0; l < 2; l++)
          coeff += beta(j,l) * beta(i,k) * PK2(k,l);
      force += p * coeff;
    }

    forceVector[2 * i] = force[0];
    forceVector[2 * i + 1] = force[1];
  }
  forceVector *= restArea();
  forceVector *= -1.0;
  return forceVector;
}
*/

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE::computeForceVectorCauchy()
{
  MATRIX2 F = computeF();
  MATRIX2 cauchy = _material->cauchy(F);
  const Real J = F.determinant();

  MATRIX2 P = J * (cauchy * F.transpose().inverse());

  MATRIX2 PK1 = _material->PK1(F);

  const Real diff = (P - PK1).norm();
  if (diff > 1e-4)
  {
    cout << " Diff: " << diff << endl;
    cout << " PK1: \n" << PK1 << endl;
    cout << " From Cauchy: \n" << P << endl;
  }

  VECTOR forceVector(6);
  VECTOR flat = flatten(P);
  forceVector = pFpu().transpose() * flat;

  return -1.0 * forceVector * restArea();
}

///////////////////////////////////////////////////////////////////////
// what are the min and max angles on the interior of this triangle?
///////////////////////////////////////////////////////////////////////
void TRIANGLE::minMaxAngles(Real& minDegrees, Real& maxDegrees)
{
  const Real degrees0 = interiorAngle(0);
  minDegrees = degrees0;
  maxDegrees = degrees0;

  for (int x = 1; x < 3; x++)
  {
    Real degrees = interiorAngle(x);

    minDegrees = (degrees < minDegrees) ? degrees : minDegrees;
    maxDegrees = (degrees > maxDegrees) ? degrees : maxDegrees;
  }
}

///////////////////////////////////////////////////////////////////////
// what is the interior angle at this vertex?
///////////////////////////////////////////////////////////////////////
Real TRIANGLE::interiorAngle(const int index) const
{
  const int next = (index + 1) % 3;
  const int nextNext = (index + 2) % 3;

  const VEC2 v0 = (_restPose[next] - _restPose[index]).normalized();
  const VEC2 v1 = (_restPose[nextNext] - _restPose[index]).normalized();

  const Real radians = acos(v0.dot(v1));

  return radians * 180.0 / M_PI;
}

///////////////////////////////////////////////////////////////////////
// build an equilateral triangle with specific edge lengths
///////////////////////////////////////////////////////////////////////
vector<VEC2> TRIANGLE::equilateral(const Real edgeLength)
{
  vector<VEC2> vertices;
  vertices.push_back(VEC2(0,0));
  vertices.push_back(VEC2(edgeLength,0));
  vertices.push_back(VEC2(edgeLength * 0.5, sqrt(3.0) * 0.5 * edgeLength));
  return vertices;
}

///////////////////////////////////////////////////////////////////////
// build the area score from [Knupp 2003]
///////////////////////////////////////////////////////////////////////
Real TRIANGLE::areaQuality()
{
  // find the longest edge length
  Real edgeLengths[3];
  for (int x = 0; x < 3; x++)
  {
    int next = (x + 1) % 3;
    VEC2 diff = _restPose[x] - _restPose[next];
    edgeLengths[x] = diff.norm();
  }
  Real maxLength = edgeLengths[0];
  for (int x = 1; x < 3; x++)
    maxLength = (maxLength > edgeLengths[x]) ? maxLength : edgeLengths[x];

  // build a canonical equilateral triangle and its Dm
  vector<VEC2> canonical = equilateral(maxLength);
  MATRIX2 T;
  T(0,0) = canonical[1][0] - canonical[0][0];
  T(1,0) = canonical[1][1] - canonical[0][1];
  T(0,1) = canonical[2][0] - canonical[0][0];
  T(1,1) = canonical[2][1] - canonical[0][1];
  Real areaCanonical = T.determinant() * 0.5;

  MATRIX2 F;
  F(0,0) = _restPose[1][0] - _restPose[0][0];
  F(1,0) = _restPose[1][1] - _restPose[0][1];
  F(0,1) = _restPose[2][0] - _restPose[0][0];
  F(1,1) = _restPose[2][1] - _restPose[0][1];
  Real area = F.determinant() * 0.5;

  // compare the areas
  Real tau = area / areaCanonical;

  // shouldn't actually need this min if we're using
  // the longest edge
  return min(tau, 1.0 / tau);
}

///////////////////////////////////////////////////////////////////////
// build the angle score from [Knupp 2003]
///////////////////////////////////////////////////////////////////////
Real TRIANGLE::angleQuality()
{
  MATRIX2 F;
  F(0,0) = _restPose[1][0] - _restPose[0][0];
  F(1,0) = _restPose[1][1] - _restPose[0][1];
  F(0,1) = _restPose[2][0] - _restPose[0][0];
  F(1,1) = _restPose[2][1] - _restPose[0][1];
  Real alpha = F.determinant();

  MATRIX2 lambda = F.transpose() * F;

  return sqrt(3.0) * alpha / (lambda(0,0) + lambda(1,1) - lambda(0,1));
}

///////////////////////////////////////////////////////////////////////
// get the maximum singular value of the material matrix
///////////////////////////////////////////////////////////////////////
Real TRIANGLE::minDmSingularValue()
{
  MATRIX2 F;
  F(0,0) = _restPose[1][0] - _restPose[0][0];
  F(1,0) = _restPose[1][1] - _restPose[0][1];
  F(0,1) = _restPose[2][0] - _restPose[0][0];
  F(1,1) = _restPose[2][1] - _restPose[0][1];

  JacobiSVD<MATRIX2> svd(F);
  return svd.singularValues()[1];
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real TRIANGLE::cTcPartial(const int index, const int goodEdge)
{
  int zero = goodEdge;
  int one = (zero + 1) % 3;
  
  const VEC2& v0 = *_vertices[zero];
  const VEC2& v1 = *_vertices[one];

  switch (index) {
    case 0:
      return 2.0 * (v0[0] - v1[0]);
    case 1:
      return 2.0 * (v0[1] - v1[1]);
    case 2:
      return 2.0 * (-v0[0] + v1[0]);
    case 3:
      return 2.0 * (-v0[1] + v1[1]);
    default:
      return 0;
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX2 TRIANGLE::CPartial(const int index, const int goodEdge)
{
  int zero = goodEdge;
  int one = (zero + 1) % 3;

  const VEC2& v0 = *_vertices[zero];
  const VEC2& v1 = *_vertices[one];

  MATRIX2 result = MATRIX2::Zero();
  switch (index) {
    case 0:
      result << 2.0 * (v0[0] - v1[0]), v0[1] - v1[1],
                       v0[1] - v1[1],           0.0;      
      break;
    case 1:
      result <<           0.0,        v0[0] - v1[0],
                v0[0] - v1[0], 2.0 * (v0[1] - v1[1]);
      break;
    case 2:
      result << 2.0 * (-v0[0] + v1[0]), -v0[1] + v1[1],
                       -v0[1] + v1[1],            0.0;      
      break;
    case 3:
      result <<            0.0,        -v0[0] + v1[0],
                -v0[0] + v1[0], 2.0 * (-v0[1] + v1[1]);
      break;
    default:
      break;
  }
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX2 TRIANGLE::DsPartial(const int index)
{
  MATRIX2 result = MATRIX2::Zero();

  if (index == 0)
  {
    result(0,0) = -1;
    result(0,1) = -1;
  }
  if (index == 1)
  {
    result(1,0) = -1;
    result(1,1) = -1;
  }
  if (index == 2)
    result(0,0) = 1;
  if (index == 3)
    result(1,0) = 1;
  if (index == 4)
    result(0,1) = 1;
  if (index == 5)
    result(1,1) = 1;

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE::pFpuDegenerateComplete(const MATRIX2& DmInverse)
{
  const int goodEdge = findLeastDegenerateEdge();
  int zero = goodEdge;
  int one = (zero + 1) % 3;
  int two = (zero + 2) % 3;

  // get rest and deformed clean direction
  VEC2 c    = *_vertices[one] - *_vertices[zero];
  VEC2 cBar = _restPose[one] - _restPose[zero];
  MATRIX2 C = c * c.transpose();
  Real cTc = c.dot(c);

  // get Dm
  VEC2 edgeBar0 = _restPose[one] - _restPose[zero];
  VEC2 edgeBar2 = _restPose[two] - _restPose[zero];
  MATRIX2 Dm;
  Dm.col(0) = edgeBar0;
  Dm.col(1) = edgeBar2;

  // get Ds
  VEC2 edge0 = *_vertices[one] - *_vertices[zero];
  VEC2 edge2 = *_vertices[two] - *_vertices[zero];
  MATRIX2 Ds;
  Ds.col(0) = edge0;
  Ds.col(1) = edge2;

  // get the constant matrices
  MATRIX2 Delta0;
  Delta0 << 1,0,0,0;
  MATRIX2 Delta1;
  Delta1 << 0,0,0,1;
  MATRIX2 E;
  E << 0,1,0,0;
  MATRIX2 T;
  T << 0,-1,1,0;

  MATRIX result(4,6);
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " Dm: " << endl << Dm << endl;
  //cout << " Ds: " << endl << Ds << endl;
  for (int x = 0; x < 6; x++)
  {
    MATRIX2 column;
    column.setZero();

    Real cTcGradient = cTcPartial(x, goodEdge);
    MATRIX2 CGradient = CPartial(x, goodEdge);
    MATRIX2 DsGradient = DsPartial(x);

    column += cTcGradient * (-1.0 / (cTc * cTc) * C * (Ds * Delta0 + Dm * Delta1));
    column += (1.0 / cTc) * CGradient * (Ds * Delta0 + Dm * Delta1);
    column += (1.0 / cTc) * C * DsGradient * Delta0;
    column += -0.5 * cTcGradient * (cBar.norm() / pow(cTc, 1.5)) * T * Ds * E;
    column += (cBar.norm() / pow(cTc, 0.5)) * T * DsGradient * E;

    //cout << " Column " << x << ": " << endl << column << endl;

    column = column * DmInverse;
    /*
    if (x == 0) 
    {
      cout << endl << column << endl;
    }
    */
    result.col(x) = flatten(column);
  }
  return result;
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

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
bool TRIANGLE::testDegenerateDFDu()
{
  MATRIX2 DmIsotropic;
  MATRIX2 DsIsotropic;
  MATRIX2 DsAnisotropic;
  VEC2 cBar;
  Real phantomRestArea;
  MATRIX6 permutation;
#if USING_ANISOTROPIC_ARAP
  //ANISOTROPIC_ARAP anisotropic(SUPER_GLUE, SUPER_GLUE);
  ANISOTROPIC_ARAP anisotropic(superGlue, superGlue);
#else
  //ANISOTROPIC_DIRICHLET anisotropic(SUPER_GLUE, SUPER_GLUE);
  ANISOTROPIC_DIRICHLET anisotropic(superGlue, superGlue);
#endif
  computeDegenerateQuantities(DmIsotropic, 
                              DsIsotropic,
                              DsAnisotropic,
                              cBar, 
                              phantomRestArea, 
                              &anisotropic,
                              permutation);
  MATRIX2 F0 = DsIsotropic * DmIsotropic.inverse();

  MATRIX pFpuDirect = pFpuDegenerateComplete(DmIsotropic.inverse());
  VECTOR uOriginal = getDisplacements();
  cout << " Direct pPpu: " << endl;
  cout << pFpuDirect << endl;

  // make a finite diff approximation
  Real eps = 1e-2;
  Real diff = 0;
  MATRIX pFpuNumerical = pFpuDirect;
  pFpuNumerical.setZero();
  for (int e = 0; e < 6; e++)
  {
    for (int i = 0; i < 6; i++)
    {
      VECTOR u = uOriginal;
      u[i] += eps;
      setDisplacements(u);
  
      computeDegenerateQuantities(DmIsotropic, 
                                  DsIsotropic,
                                  DsAnisotropic,
                                  cBar, 
                                  phantomRestArea, 
                                  &anisotropic,
                                  permutation);
      MATRIX2 F1 = DsIsotropic * DmIsotropic.inverse();
      MATRIX2 finiteDiff = F1 - F0;
      finiteDiff *= 1.0 / eps;

      pFpuNumerical.col(i) = flatten(finiteDiff);
    }
    //MATRIX pFpuDiff = DPDF - Pnumerical;
    MATRIX pFpuDiff = pFpuDirect - pFpuNumerical;
    diff = pFpuDiff.norm() / pFpuDiff.size();

    cout << " eps: \t" << eps << " \t diff: \t" << diff << endl;
#if 1
    if (e == 5)
    {
      cout << " numerical: " << endl;
      cout << pFpuNumerical << endl;
      cout << " diff: " << endl;
      cout << pFpuDiff << endl;
    }
#endif

    eps *= 0.1;
  }
  cout << " Non-degenerate pFpu: " << endl;
  cout << pFpu() << endl;

  return true;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE::runDegeneracyUnitTests()
{
  // generate the vertices
  vector<VEC2> vertices;
  vertices.push_back(VEC2(0.0, 0.0));
  vertices.push_back(VEC2(0.0, 1.0));
  //vertices.push_back(VEC2(0.0, 0.5678));
  vertices.push_back(VEC2(-1e-4, 0.5678));
  //vertices.push_back(VEC2(-1e-9, 0.5));

  // get the pointers
  vector<VEC2*> vertexPointers;
  for (int x = 0; x < 3; x++)
    vertexPointers.push_back(&vertices[x]);

  // build the triangle
  ARAP arap(1.0, 1.0);
  TRIANGLE triangle(&arap, vertexPointers);

  vertices[0][0] -= 1.0;
  vertices[1][0] -= 1.0;
  vertices[2][0] += 1.0;

  triangle.testDegenerateDFDu();
}
