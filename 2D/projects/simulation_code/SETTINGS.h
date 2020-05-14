#ifndef SETTINGS_H
#define SETTINGS_H

#include <cassert>

#include "Eigen/Dense"
using namespace Eigen;

#define Real double

typedef Matrix<Real, 2, 2 > MATRIX2;
typedef Matrix<Real, 3, 3 > MATRIX3;
typedef Matrix<Real, 4, 4 > MATRIX4;
typedef Matrix<Real, 6, 6 > MATRIX6;
typedef Matrix<Real, 2, 1 > VEC2;
typedef Matrix<Real, 3, 1 > VEC3;
typedef Matrix<Real, 4, 1 > VEC4;
typedef Matrix<int, 3, 1 > VEC3I;

typedef Matrix<Real, Dynamic, Dynamic> MATRIX;
typedef Matrix<Real, Dynamic, 1> VECTOR;

#endif
