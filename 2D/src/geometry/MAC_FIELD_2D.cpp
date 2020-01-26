#include "MAC_FIELD_2D.h"
#include <GL/freeglut.h>

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MAC_FIELD_2D::MAC_FIELD_2D(const int& xRes, const int& yRes, const VEC3& center, const VEC3& lengths) :
  _xRes(xRes), _yRes(yRes), _center(center), _lengths(lengths),
  _pressure(xRes, yRes, center, lengths)
{
  _totalCells = _xRes * _yRes;
  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;

  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;

  // allocate the x and y velocity fields
  _xVelocity = FIELD_2D(xRes + 1, _yRes, _center, _lengths + VEC3(_dx,0,0));
  _yVelocity = FIELD_2D(xRes, _yRes + 1, _center, _lengths + VEC3(0,_dy,0));
}

MAC_FIELD_2D::MAC_FIELD_2D() :
  _xRes(-1), _yRes(-1), _totalCells(-1), _dx(-1), _dy(-1),
  _invDx(-1), _invDy(-1)
{
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MAC_FIELD_2D::~MAC_FIELD_2D()
{
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void MAC_FIELD_2D::drawGrid() const
{
  VEC3 lowerLeft  = _center - 0.5 * _lengths; 
  VEC3 upperRight = _center + 0.5 * _lengths; 

  // draw the bounding box
  glColor4f(1,1,1,1);

  glBegin(GL_LINE_STRIP);
    glVertex2f(lowerLeft[0], lowerLeft[1]);
    glVertex2f(lowerLeft[0], upperRight[1]);
    glVertex2f(upperRight[0], upperRight[1]);
    glVertex2f(upperRight[0], lowerLeft[1]);
    glVertex2f(lowerLeft[0], lowerLeft[1]);
  glEnd();

  // draw the horizontal grid lines
  //glColor4f(0,1,0,1);
  Real yCurrent = lowerLeft[1] + _dy;
  glBegin(GL_LINES);
  for (int x = 0; x < _xRes - 1; x++)
  {
    glVertex2f(lowerLeft[0], yCurrent);
    glVertex2f(upperRight[0], yCurrent);
    yCurrent += _dy;
  }
  glEnd();

  // draw the vertical grid lines
  //glColor4f(0,0,1,1);
  Real xCurrent = lowerLeft[0] + _dx;
  glBegin(GL_LINES);
  for (int x = 0; x < _xRes - 1; x++)
  {
    glVertex2f(xCurrent, lowerLeft[1]);
    glVertex2f(xCurrent, upperRight[1]);
    xCurrent += _dx;
  }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void MAC_FIELD_2D::drawMarkers() const
{
  // draw the pressure centers
  glColor4f(0,1,0,1);
  _pressure.drawCellCenters();
  
  // draw the velocity centers
  glColor4f(1,0,0,1);
  _xVelocity.drawCellCenters();
  glColor4f(0,0,1,1);
  _yVelocity.drawCellCenters();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void MAC_FIELD_2D::drawCurlVectors() const
{
  const int xRes = _pressure.xRes();
  const int yRes = _pressure.yRes();

  glBegin(GL_LINES);
  for (int y = 0; y < yRes; y++)
    for (int x = 0; x < xRes; x++)
    {
      VEC3 center = _pressure.cellCenter(x,y);

      Real xVelocity = _xVelocity(center);
      Real yVelocity = _yVelocity(center);

      Vector3d v(0,0,1);
      Vector3d w(xVelocity, yVelocity, 0);

      Vector3d curl = v.cross(w);

      glColor4f(1,1,0,1);
      glVertex2f(center[0], center[1]);
      glColor4f(1,1,0,0);
      //glVertex2f(center[0] + xVelocity, center[1] + yVelocity);
      glVertex2f(center[0] + curl[0], center[1] + curl[1]);


      /*
      if (y == 0)
      {
        cout << " center: " << center << endl;
        cout << " xy: " << x << " " << y << " velocity: " << xVelocity << " " << yVelocity << endl;
      } 
     */ 
    }
  glEnd();
  //exit(0);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void MAC_FIELD_2D::drawVectors() const
{
  const int xRes = _pressure.xRes();
  const int yRes = _pressure.yRes();

  glBegin(GL_LINES);
  for (int y = 0; y < yRes; y++)
    for (int x = 0; x < xRes; x++)
    {
      VEC3 center = _pressure.cellCenter(x,y);

      Real xVelocity = _xVelocity(center);
      Real yVelocity = _yVelocity(center);

      glColor4f(1,1,0,1);
      glVertex2f(center[0], center[1]);
      glColor4f(1,1,0,0);
      glVertex2f(center[0] + xVelocity, center[1] + yVelocity);


      /*
      if (y == 0)
      {
        cout << " center: " << center << endl;
        cout << " xy: " << x << " " << y << " velocity: " << xVelocity << " " << yVelocity << endl;
      } 
     */ 
    }
  glEnd();
  //exit(0);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D MAC_FIELD_2D::divergence()
{
  // final result has the same dimensions as the pressure field
  FIELD_2D final(_pressure);

  int yRes = _pressure.yRes();
  int xRes = _pressure.xRes();
  for (int y = 0; y < yRes; y++)
    for (int x = 0; x < xRes; x++)
      final(x,y) = (_xVelocity(x,y) - _xVelocity(x + 1, y)) * _invDx +
                   (_yVelocity(x,y) - _yVelocity(x, y + 1)) * _invDy;

  return final;
}

///////////////////////////////////////////////////////////////////////
// set the velocity field to a swirl
///////////////////////////////////////////////////////////////////////
void MAC_FIELD_2D::setToSwirl()
{
  for (int y = 0; y < _xVelocity.yRes(); y++)
    for (int x = 0; x < _xVelocity.xRes(); x++)
    {
      VEC3 center = _xVelocity.cellCenter(x,y);
      center -= _center;

      center[0] *= 1.0 / _lengths[0];
      center[1] *= 1.0 / _lengths[1];
      center *= 2.0 * M_PI;

      //_xVelocity(x,y) =  sin(center[0]) * cos(center[1]);
      //_yVelocity(x,y) = -cos(center[0]) * sin(center[1]);
      _xVelocity(x,y) =  sin(center[0]) * cos(0 * center[1]);
      _yVelocity(x,y) = -cos(0 * center[0]) * sin(center[1]);

      /*
      if (y == 0)
      {
        cout << " xy: " << x << " " << y << endl;
        cout << " center: " << center << endl;
        cout << " sine: " << sin(center[0]) << endl;
      }
      */
    }

  _xVelocity *= 0.1;
  _yVelocity *= 0.1;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void MAC_FIELD_2D::advect(const Real& dt, const FIELD_2D& oldField, FIELD_2D& newField)
{
  assert(oldField.xRes() == newField.xRes());  
  assert(oldField.yRes() == newField.yRes());  

  const int yRes = oldField.yRes();
  const int xRes = oldField.xRes();
  for (int y = 0; y < yRes; y++)
    for (int x = 0; x < xRes; x++)
    {
      VEC3 center = oldField.cellCenter(x,y);
      VEC3 velocity(_xVelocity(center), _yVelocity(center), 0);

      // trace the ray backwards
      VEC3 backtrace = center - dt * velocity;
      newField(x,y) = oldField(backtrace);
    }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void MAC_FIELD_2D::advectMacCormack(const Real& dt, const FIELD_2D& oldField, FIELD_2D& newField)
{
	FIELD_2D phiHatN(oldField);
	FIELD_2D phiHatN1(oldField);

	for (int x = 0; x < oldField.totalCells(); x++)
		phiHatN[x] = phiHatN1[x] = oldField[x];

	const FIELD_2D& phiN = oldField;
	FIELD_2D& phiN1 = newField;

	// phiHatN1 = A(phiN)
  cout << " forward ... ";flush(cout);
	advect(dt, phiN, phiHatN1);

	// phiHatN = A^R(phiHatN1)
  cout << " backward ... ";flush(cout);
	advect(-1.0 * dt, phiHatN1, phiHatN);

  phiN1.clear();
  phiN1 += phiN;
  phiN1 -= phiHatN;
  phiN1 *= 0.5;
  phiN1 += phiHatN1;

	// clamp any newly created extrema
  cout << " clamping extrema ... ";flush(cout);
	clampExtrema(dt, oldField, newField);

	// if the error estimate was bad, revert to first order
  cout << " clamping rays ... ";flush(cout);
	clampOutsideRays(dt, oldField, phiHatN1, newField);
  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
// Clamp the extrema generated by the BFECC error correction
///////////////////////////////////////////////////////////////////////
void MAC_FIELD_2D::clampExtrema(const Real& dt, const FIELD_2D& oldField, FIELD_2D& newField)
{
  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();

  for (int y = 1; y < yRes-1; y++)
		for (int x = 1; x < xRes-1; x++)
    {
      // backtrace
      const VEC3 center = oldField.cellCenter(x,y);
      const VEC3 velocity(_xVelocity(center), _yVelocity(center), 0);
      if (velocity.norm() < 1e-6) continue;

      Real xTrace = x - dt * velocity[0];
      Real yTrace = y - dt * velocity[1];
      //Real zTrace = z - dt * velocity[2];

      // clamp backtrace to grid boundaries
      xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
      xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
      yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
      yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
      //zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
      //zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

      // locate neighbors to interpolate
      const int x0 = (int)xTrace;
      const int x1 = x0 + 1;
      const int y0 = (int)yTrace;
      const int y1 = y0 + 1;
      //const int z0 = (int)zTrace;
      //const int z1 = z0 + 1;

      Real minField = oldField(x0,y0);
      Real maxField = minField;

      const Real x0y1z0 = oldField(x0,y1);
      const Real x1y0z0 = oldField(x1,y0);
      const Real x1y1z0 = oldField(x1,y1);

      minField = (x0y1z0 < minField) ? x0y1z0 : minField;
      maxField = (x0y1z0 > maxField) ? x0y1z0 : maxField;

      minField = (x1y0z0 < minField) ? x1y0z0 : minField;
      maxField = (x1y0z0 > maxField) ? x1y0z0 : maxField;

      minField = (x1y1z0 < minField) ? x1y1z0 : minField;
      maxField = (x1y1z0 > maxField) ? x1y1z0 : maxField;

      const Real newValue = newField(x,y);

      newField(x,y) = (newValue > maxField) ? maxField : newValue;
      newField(x,y) = (newValue < minField) ? minField : newValue;
    }
}

///////////////////////////////////////////////////////////////////////
// Reverts any backtraces that go into boundaries back to first 
// order -- in this case the error correction term was totally
// incorrect
///////////////////////////////////////////////////////////////////////
void MAC_FIELD_2D::clampOutsideRays(const Real& dt, const FIELD_2D& oldField, const FIELD_2D& oldAdvection, FIELD_2D& newField)
{
  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();

  for (int y = 0; y < yRes; y++)
		for (int x = 0; x < xRes; x++)
    {
      const VEC3 center = oldField.cellCenter(x,y);
      VEC3 velocity(_xVelocity(center), _yVelocity(center), 0);
      if (velocity.norm() < 1e-6) continue;
      velocity[0] *= dt;
      velocity[1] *= dt;

      float xBackward = x + velocity[0];
      float yBackward = y + velocity[1];

      float xTrace    = x - velocity[0];
      float yTrace    = x - velocity[1];

      // see if it goes outside the boundaries
      bool hasObstacle = 
        (yTrace < 1.0f)    || (yTrace > yRes - 2.0f) ||
        (xTrace < 1.0f)    || (xTrace > xRes - 2.0f) ||
        (yBackward < 1.0f) || (yBackward > yRes - 2.0f) ||
        (xBackward < 1.0f) || (xBackward > xRes - 2.0f);

      // reuse old advection instead of doing another one...
      if(hasObstacle) { newField(x,y) = oldAdvection(x,y); continue; }

    } // xyz
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void MAC_FIELD_2D::buildGradientField(const FIELD_2D& scalar)
{
  // build the x velocity field
  _xVelocity = 0;
  for (int y = 0; y < scalar.yRes(); y++)
    for (int x = 0; x < scalar.xRes() - 1; x++)
      _xVelocity(x + 1,y) = (scalar(x,y) - scalar(x+1,y)) * _invDx;
  
  // build the y velocity field
  _yVelocity = 0;
  for (int y = 0; y < scalar.yRes() - 1; y++)
    for (int x = 0; x < scalar.xRes(); x++)
      _yVelocity(x,y + 1) = (scalar(x,y) - scalar(x,y+1)) * _invDy;

  _xVelocity.unitize();
  _yVelocity.unitize();

  //_yVelocity *= -1;
  //_xVelocity *= -1;
  //_xVelocity *= 2;
  //_yVelocity *= 2;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void MAC_FIELD_2D::setToEigenvector()
{
  int k1 = 1;
  int k2 = 1;
  Real invKs = 1.0 / (k1 * k1 + k2 * k2);

  for (int y = 0; y < _yVelocity.yRes(); y++)
    for (int x = 0; x < _yVelocity.xRes(); x++)
    {
      VEC3 center = _yVelocity.cellCenter(x,y);
      _yVelocity(x,y) = -invKs * k1 * cos(k1 * center[0]) * sin(k2 * center[1]);
    }
  
  for (int y = 0; y < _xVelocity.yRes(); y++)
    for (int x = 0; x < _xVelocity.xRes(); x++)
    {
      VEC3 center = _xVelocity.cellCenter(x,y);
      _xVelocity(x,y) = invKs * k2 * sin(k1 * center[0]) * cos(k2 * center[1]);
    }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VectorXd MAC_FIELD_2D::flattenVelocity() const
{
  int xCells = _xVelocity.totalCells();
  int yCells = _yVelocity.totalCells();
  int cells = xCells + yCells;
  VectorXd final(cells);

  for (int x = 0; x < xCells; x++)
    final[x] = _xVelocity[x];
  for (int y = 0; y < yCells; y++)
    final[y + xCells] = _yVelocity[y];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void MAC_FIELD_2D::unflattenVelocity(const VectorXd& input)
{
  int xCells = _xVelocity.totalCells();
  int yCells = _yVelocity.totalCells();

  for (int x = 0; x < xCells; x++)
    _xVelocity[x] = input[x];
  for (int y = 0; y < yCells; y++)
    _yVelocity[y] = input[y + xCells];
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void MAC_FIELD_2D::unflattenVelocity(const VectorXf& input)
{
  int xCells = _xVelocity.totalCells();
  int yCells = _yVelocity.totalCells();

  for (int x = 0; x < xCells; x++)
    _xVelocity[x] = input[x];
  for (int y = 0; y < yCells; y++)
    _yVelocity[y] = input[y + xCells];
}
