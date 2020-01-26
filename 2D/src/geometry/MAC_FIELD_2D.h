#ifndef MAC_FIELD_2D_H
#define MAC_FIELD_2D_H

#include <cmath>
#include <string>
#include <iostream>
#include <SETTINGS.h>

#include "FIELD_2D.h"

using namespace std;

class MAC_FIELD_2D {
public:
  MAC_FIELD_2D();
  MAC_FIELD_2D(const int& xRes, const int& yRes, const VEC3& center = VEC3(0,0,0), const VEC3& lengths = VEC3(1,1,1));
  MAC_FIELD_2D(const MAC_FIELD_2D& m);
  ~MAC_FIELD_2D();

  void drawGrid() const;
  void drawMarkers() const;
  void drawVectors() const;
  void drawCurlVectors() const;

  // accessors
  /*
  inline Real& operator()(int x, int y) { return _data[y * _xRes + x]; };
  inline Real operator()(int x, int y) const { return _data[y * _xRes + x]; };
  inline Real& operator[](int x) { return _data[x]; };
  inline Real operator[](int x) const { return _data[x]; };
  inline Real& operator()(const VEC3& v) { return _data[cellIndex(v)]; };
  Real* data() { return _data; };
  */
  int xRes() const { return _xRes; };
  int yRes() const { return _yRes; };
  const VEC3& center() const { return _center; };
  const VEC3& lengths() const { return _lengths; };
  int totalCells() const { return _totalCells; };
  const Real& dx() const { return _dx; };
  const Real& dy() const { return _dy; };
  const Real& invDx() const { return _invDx; };
  const Real& invDy() const { return _invDy; };
  VEC3 dims() const { return VEC3(_xRes, _yRes, 0); };

  FIELD_2D& pressure() { return _pressure; };
  FIELD_2D& xVelocity() { return _xVelocity; };
  FIELD_2D& yVelocity() { return _yVelocity; };

  FIELD_2D divergence();

  void advect(const Real& dt, const FIELD_2D& oldField, FIELD_2D& newField);
  void advectMacCormack(const Real& dt, const FIELD_2D& oldField, FIELD_2D& newField);
  void clampExtrema(const Real& dt, const FIELD_2D& oldField, FIELD_2D& newField);
  void clampOutsideRays(const Real& dt, const FIELD_2D& oldField, const FIELD_2D& oldAdvection, FIELD_2D& newField);

  void setToSwirl();

  void buildGradientField(const FIELD_2D& scalar);

  void setToEigenvector();

  VectorXd flattenVelocity() const;
  void unflattenVelocity(const VectorXd& input);
  void unflattenVelocity(const VectorXf& input);

  /*
  // get the real-space location of a cell
  VEC3 cellCenter(int index) const;
  VEC3 cellCenter(int x, int y) const;
  int cellIndex(const VEC3& position) const;
  void cellIndex(const VEC3& position, int& x, int& y) const;

  // a safe, slow, toroidal data accessor -- does all bounds checking for you
  inline Real safe(int x, int y);

  // common field operations
  void clear();
  void normalize();
  MAC_FIELD_2D& abs();

  Real min();
  Real max();

  // field maximum cell index
  VEC3 maxIndex();

  // field minimum cell index
  VEC3 minIndex();

  // take the log
  void log(Real base = 2.0);
 
  // IO functions
  void writeJPG(string filename);
  void writePPM(string filename);
  void writeMatlab(string filename, string variableName) const;
  void write(string filename) const;
  void write(FILE* file) const;
  void read(string filename);
  void read(FILE* file);
  void writePNG(string filename);
  void readPNG(string filename);

  // FFT functions
  void FFT(MAC_FIELD_2D& real, MAC_FIELD_2D& im);
  void shiftFFT();

  // set this field to the result of convolving filter and input
  void convolve(const MAC_FIELD_2D& filter, const MAC_FIELD_2D& input);

  void resizeAndWipe(int xRes, int yRes, const VEC3& center, const VEC3& lengths);

  // overloaded operators
  MAC_FIELD_2D& operator=(const Real& alpha);
  MAC_FIELD_2D& operator=(const MAC_FIELD_2D& A);
  MAC_FIELD_2D& operator*=(const Real& alpha);
  MAC_FIELD_2D& operator/=(const Real& alpha);
  MAC_FIELD_2D& operator+=(const Real& alpha);
  MAC_FIELD_2D& operator-=(const Real& alpha);
  MAC_FIELD_2D& operator-=(const MAC_FIELD_2D& input);
  MAC_FIELD_2D& operator+=(const MAC_FIELD_2D& input);
  MAC_FIELD_2D& operator*=(const MAC_FIELD_2D& input);
  MAC_FIELD_2D& operator/=(const MAC_FIELD_2D& input);

  // sum of all entries
  Real sum();
  
  // Compute the elements of the vertical derivative convolution kernel
  void verticalDerivativeKernel(double kMax = 10, double dk = 0.01, double sigma = 1.0, double L = 0);
  
  // Compute a radial Bessel function
  void radialBessel();
  
  // set to a bessel function
  void setToBessel(float k);

  MAC_FIELD_2D nearestNeighborUpsample(int factor);

  // set to a checkboard for debugging
  void setToCheckerboard(int xChecks = 10, int yChecks = 10);
  
  // set to a checkboard for debugging
  void setToRampedCheckerboard(int xChecks = 10, int yChecks = 10);
  
  // set to a ramp for debugging
  void setToRampX();
  
  // set to a ramp for debugging
  void setToRampY();

  // pass a field to fieldViewer2D
  static void fieldViewer(const MAC_FIELD_2D& field, string name);
  static void overlayFieldViewer(const MAC_FIELD_2D& distanceField, const MAC_FIELD_2D& field, string name);

  // get the projection of the field in the x direction
  //VECTOR projectionX();

  // return a field for the Laplacian of this field
  MAC_FIELD_2D laplacian();
  
  // return a field for the 4th order accurate Laplacian of this field
  MAC_FIELD_2D laplacian4th();
  
  // return a field for the gradient of this field
  MAC_FIELD_2D gradient();

  // return the transpose (flip x and y)
  MAC_FIELD_2D transpose() const;

  // use wavelet upsampling to double resolution
  MAC_FIELD_2D doubleRes() const;

  // get the sum of the absolute value of all the entries
  Real fabsSum();

  // get the mean curvature
  MAC_FIELD_2D meanCurvature() const;
  MAC_FIELD_2D gaussianCurvature() const;

  // first order spatial derivatives
  // on the border, difference falls back to first order (not centered) difference
  Real Dx(int x, int y) const; 
  Real Dy(int x, int y) const; 

  // second order spatial derivatives
  // on the border, center cell is copied to outside, and centered difference
  // is still taken
  Real DDx(int x, int y) const; 
  Real DDy(int x, int y) const; 

  // mixed derivatives
  // on the border, center cell is copied to outside, and centered difference
  // is still taken
  Real DDxy(int x, int y) const;

  // get the normal at a point
  VEC3 normal(int x, int y) const;

  // count the number of infs in the field
  int totalInfs() const;
  
  // count the number of NaNs in the field
  int totalNans() const;
  */

private:
  int _xRes;
  int _yRes;
  int _totalCells;

  VEC3 _center;
  VEC3 _lengths;
  Real _dx, _dy;
  Real _invDx, _invDy;

  FIELD_2D _xVelocity;
  FIELD_2D _yVelocity;
  FIELD_2D _pressure;
};

/*
MAC_FIELD_2D operator*(const MAC_FIELD_2D& A, const Real alpha);
MAC_FIELD_2D operator/(const MAC_FIELD_2D& A, const Real alpha);
MAC_FIELD_2D operator+(const MAC_FIELD_2D& A, const Real alpha);
MAC_FIELD_2D operator*(const Real alpha, const MAC_FIELD_2D& A);
MAC_FIELD_2D operator+(const Real alpha, const MAC_FIELD_2D& A);
MAC_FIELD_2D operator-(const MAC_FIELD_2D& A, const MAC_FIELD_2D& B);
MAC_FIELD_2D operator+(const MAC_FIELD_2D& A, const MAC_FIELD_2D& B);
MAC_FIELD_2D operator*(const MAC_FIELD_2D& A, const MAC_FIELD_2D& B);
*/

#endif
