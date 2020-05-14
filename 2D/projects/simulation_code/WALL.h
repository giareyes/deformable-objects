#ifndef WALL_H
#define WALL_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "SETTINGS.h"
using namespace std;

class WALL {
public:
  WALL(const VEC2& normal, const VEC2& point);

  // accessors
  VEC2& normal() { return _normal; };
  VEC2& point()  { return _point; };

private:
  VEC2 _normal;
  VEC2 _point;
};

#endif
