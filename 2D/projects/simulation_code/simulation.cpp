#include "SETTINGS.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "TRIANGLE_MESH.h"
#include "STVK.h"
#include "NEOHOOKEAN.h"
#include "WALL.h"

using namespace std;

// the resolution of the OpenGL window -- independent of the field resolution
int xScreenRes = 800;
int yScreenRes = 800;

// Text for the title bar of the window
string windowLabel("Reduced-Order Neo-Hookean Sim");

// mouse tracking variables
int xMouse         = -1;
int yMouse         = -1;
int mouseButton    = -1;
int mouseState     = -1;
int mouseModifiers = -1;

// animate the current runEverytime()?
bool animate = false;
bool singleStep = false;
bool createBasis = false;
// float dt = 1.0/15360.0;
float dt = 1.0/100.0;

// the current viewer eye position
VEC3 eyeCenter(0, 0, 1);

// current zoom level into the field
float zoom = 2.0;

//Real poissonsRatio = 0.0;
// old:
// Real poissonsRatio = 0.4;
// Real youngsModulus = 1.0;
int basisCols = 0;

// testing:
Real poissonsRatio = 0.49;
Real youngsModulus = 1.0;

TRIANGLE_MESH triangleMesh(poissonsRatio, youngsModulus);
VEC2 bodyForce;

enum SCENE { STRETCH, SQUASH, LSHEAR, RSHEAR, SINGLE, MOTION};
SCENE scene = SINGLE;
int sceneNum = 1;

int meshFlag = 0;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawTriangle(const TRIANGLE& triangle, const VEC4& color)
{
  glColor4f(color[0], color[1], color[2], color[3]);
  glBegin(GL_TRIANGLES);
    for (int x = 0; x < 3; x++)
      glVertex2f(triangle.vertex(x)[0], triangle.vertex(x)[1]);
  glEnd();

  if(meshFlag == 1)
  {
    glColor4f(0,0,0,1);
    glBegin(GL_LINE_LOOP);
      for (int x = 0; x < 3; x++)
        glVertex2f(triangle.vertex(x)[0], triangle.vertex(x)[1]);
  }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawMesh(const VEC4& color1, const VEC4& color2)
{
  // draw the triangles
  const vector<TRIANGLE>& triangles = triangleMesh.triangles();
  for (unsigned int x = 0; x < triangles.size(); x++)
  {
    drawTriangle(triangles[x], color1);
  }


  // get vertices
  vector<VEC2>& vertices = triangleMesh.vertices();
  const vector<int>& constrainedVertices = triangleMesh.constrainedVertices();

  // draw the constrained vertices if mesh flag is on
  if(meshFlag == 1)
  {
    glPointSize(5.0);
    glColor4f(1,0,0,1);
    glBegin(GL_POINTS);
    for (unsigned int x = 0; x < constrainedVertices.size(); x++)
    {
      glVertex2f(vertices[constrainedVertices[x]][0],
                 vertices[constrainedVertices[x]][1]);
    }
  }
  glEnd();
}

void drawWalls()
{
  vector<WALL>& walls = triangleMesh.walls();
  // draw the walls
  glColor4f(101.0/255,106.0/255,110.0/255,1);
  for (unsigned int x = 0; x < walls.size(); x++)
  {
    glPushMatrix();
      // translate to the point
      glTranslatef(walls[x].point()[0], walls[x].point()[1], 0);

      // apply a rotation
      float angle = asin(walls[x].normal()[0]) / (2 * M_PI) * 360.0;
      glRotatef(-angle, 0, 0, 1);

      // make it a plane at 0,0
      glTranslatef(0, -0.5, 0);
      glScalef(50,1,1);
      glutSolidCube(1.0);
    glPopMatrix();
  }
  glFlush();
}
///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  // Make ensuing transforms affect the projection matrix
  glMatrixMode(GL_PROJECTION);

  // set the projection matrix to an orthographic view
  glLoadIdentity();
  float halfZoom = zoom * 0.5;

  glOrtho(-halfZoom, halfZoom, -halfZoom, halfZoom, -10, 10);

  // set the matrix mode back to modelview
  glMatrixMode(GL_MODELVIEW);

  // set the lookat transform
  glLoadIdentity();
  gluLookAt(eyeCenter[0], eyeCenter[1], 1,  // eye
            eyeCenter[0], eyeCenter[1], 0,  // center
            0, 1, 0);   // up

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  drawMesh(VEC4(213.0 / 255.0, 246.0 / 255, 247.0 / 255.0,1.0), VEC4(157.0 / 255.0, 230.0 / 255, 156.0 / 255.0,1.0));
  drawWalls();

  glutSwapBuffers();
}

///////////////////////////////////////////////////////////////////////
// Map the keyboard keys to something here
///////////////////////////////////////////////////////////////////////
void glutKeyboard(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 'q':
      exit(0);
      break;

    case 'a':
      animate = !animate;
      break;

    case ' ':
      animate = true;
      singleStep = true;
      break;

    case 'v':
      cout << " eye: " << eyeCenter << endl;
      cout << " zoom: " << zoom << endl;
      break;

    case 's':
      break;

    default:
      break;
  }
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is clicked
///////////////////////////////////////////////////////////////////////
void glutMouseClick(int button, int state, int x, int y)
{
  int modifiers = glutGetModifiers();
  mouseButton = button;
  mouseState = state;
  mouseModifiers = modifiers;

  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && modifiers & GLUT_ACTIVE_SHIFT)
  {
    // make sure nothing else is called
    return;
  }

  xMouse = x;
  yMouse = y;
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is clicked and moving
///////////////////////////////////////////////////////////////////////
void glutMouseMotion(int x, int y)
{
  if (mouseButton == GLUT_LEFT_BUTTON &&
      mouseState == GLUT_DOWN &&
      mouseModifiers & GLUT_ACTIVE_SHIFT)
  {
    // make sure nothing else is called
    return;
  }

  vector<VEC2>& vertices = triangleMesh.vertices();
  vector<int>& unconstrainedVertices = triangleMesh.unconstrainedVertices();

  float xDiff = x - xMouse;
  float yDiff = y - yMouse;

  if (mouseButton == GLUT_LEFT_BUTTON)
  {
    // convert screen coordinates to simulation coordinates
    double screenX = ((double) x/ 400.0) - 1;
    double screenY = ((double) (800 - y) / 400.0) - 1;
    double mouseX = ((double) xMouse / 400.0) - 1;
    double mouseY = ((double) (800 - yMouse) / 400.0) - 1;

    if (sceneNum == 1)
    {
      for (int i = 0; i < unconstrainedVertices.size(); i++)
      {
        if (screenX <= vertices[unconstrainedVertices[i]][0] + 0.04
          && screenX >= vertices[unconstrainedVertices[i]][0] - 0.04
          && screenY <= vertices[unconstrainedVertices[i]][1] + 0.03
          && screenY >= vertices[unconstrainedVertices[i]][1] - 0.03)
        {
          printf("addforce\n");
          printf("vertex %d: (%f,%f)\n", i, vertices[unconstrainedVertices[i]][0], vertices[unconstrainedVertices[i]][1]);
          printf("screenx %f\n", screenX);
          printf("screeny %f\n\n\n", screenY);
          triangleMesh.addSingleForce(VEC2(xDiff, -1*yDiff), unconstrainedVertices[i]);
        }
      }
    }
    else
    {
      printf("else\n");
      for (int i = 0; i < vertices.size(); i++)
      {
        if (screenX <= vertices[i][0] + 0.04
          && screenX >= vertices[i][0] - 0.04
          && screenY <= vertices[i][1] + 0.03
          && screenY >= vertices[i][1] - 0.03)
        {
          printf("addforce\n");
          triangleMesh.addBodyForce(VEC2(xDiff*2, -1*yDiff*2));
          triangleMesh.addSingleForce(VEC2(xDiff, -1*yDiff), i);
        }
      }
    }
  }

  xMouse = x;
  yMouse = y;
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is not clicked and moving
///////////////////////////////////////////////////////////////////////
void glutPassiveMouseMotion(int x, int y)
{
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  if (animate)
  {
    static int frame = 0;
    switch (scene) {
      case STRETCH:
        triangleMesh.stretch2(0.002);
        break;
      case RSHEAR:
        triangleMesh.stepShearTest(0.002);
        break;
      case LSHEAR:
        triangleMesh.stepShearTest(-0.002);
        break;
      case SQUASH:
        triangleMesh.stretch2(-0.002);
        break;
      case SINGLE:
        triangleMesh.addBodyForce(bodyForce);
        break;
      case MOTION:
        break;
    }

    if (scene == MOTION)
    {
      if (frame == 0 && sceneNum == 0) // if we are on the first frame, insert a force so they jump
      {
        triangleMesh.addBodyForce(VEC2(10.0, 200.0));
      }
      for( int i = 0; i < 15; i ++)
      {
        if (sceneNum == 0)
        {
          triangleMesh.addBodyForce(bodyForce);
        }
        triangleMesh.stepMotion(dt, bodyForce, sceneNum);
      }
    }
    else
    {
      triangleMesh.stepQuasistatic();
    }
    frame++;

    if (singleStep)
    {
      animate = false;
      singleStep = false;
    }
  }
  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int glutWindow()
{
  glutInitDisplayMode(GLUT_DOUBLE| GLUT_RGBA);
  glutInitWindowSize(xScreenRes, yScreenRes);
  glutInitWindowPosition(10, 10);
  glutCreateWindow(windowLabel.c_str());

  // set the viewport resolution (w x h)
  glViewport(0, 0, (GLsizei) xScreenRes, (GLsizei) yScreenRes);

  // set the background color to gray
  //glClearColor(0.1, 0.1, 0.1, 0);
  glClearColor(1,1,1,1);

  // register all the callbacks
  glutDisplayFunc(&glutDisplay);
  glutIdleFunc(&glutIdle);
  glutKeyboardFunc(&glutKeyboard);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);
  glutPassiveMotionFunc(&glutPassiveMouseMotion);

  //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  //glEnable(GL_MULTISAMPLE);
  glLineWidth(1.0);
  glEnable(GL_LINE_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

  // enter the infinite GL loop
  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

string toUpper(const string& input)
{
  string copy(input);
  for (unsigned int x = 0; x < input.length(); x++)
    copy[x] = std::toupper(copy[x]);
  return copy;
}

///////////////////////////////////////////////////////////////////////
// process the command line
///////////////////////////////////////////////////////////////////////
void readCommandLine(int argc, char** argv)
{
  if (argc < 2)
  {
    cout << " invalid number of arguments. Args given: " << argc << ". need filenames" << endl;
    exit(0);
  }

  if (argc > 2)
  {
    string sceneType(argv[2]);
    sceneType = toUpper(sceneType);

    if (sceneType.compare("LSHEAR") == 0)
      scene = LSHEAR;
    else if (sceneType.compare("RSHEAR") == 0 || sceneType.compare("SHEAR") == 0 )
      scene = RSHEAR;
    else if (sceneType.compare("SQUASH") == 0)
      scene = SQUASH;
    else if (sceneType.compare("STRETCH") == 0)
      scene = STRETCH;
    else if (sceneType.compare("SINGLE") == 0)
      scene = SINGLE;
    else
    {
      scene = MOTION;
      sceneNum = 0;
    }

    if (argc > 3)
    {
      for(int x = 3; x < argc; x++)
      {
        if (argv[x][1] == 'm')
          meshFlag = 1;
        else if (argv[x][1] == 'b')
          sceneNum = 1;
        else if (argv[x][1] == 'q')
        {
          createBasis = true;
          if( x + 1 == argc )
          {
            cout << " invalid number of arguments. Need an int after the -q flag" << endl;
            exit(0);
          }
          basisCols = stoi(argv[x+1]);
          x++;
        }
      }
    }
  }
  else
  {
    scene = MOTION;
    sceneNum = 0;
  }

  // if we are creating a basis, we must create a file with all of the deformations we will be using
  // im not sure if this if=else statement will be kept later - i think ideally we should always be reducing and
  // not necessarily creating the basis file ? but this is necessary for now if we want to run quasistatics
  if(createBasis)
  {
    FILE* file = NULL;

    string filename = argv[1] + string(".basis");
    file = fopen(filename.c_str(), "w");

    for(int i = 0; i < 4; i++)
    {
      TRIANGLE_MESH basisBuild(poissonsRatio, youngsModulus);
      basisBuild.buildBlob(1, argv[1], true, 0);
      bodyForce[0] = 0;
      bodyForce[1] = -0.3;
      for(int j = 0; j < 15; j++)
      {
        switch(i) {
          case 0:
            basisBuild.stretch2(0.002);
            break;
          case 1:
            basisBuild.stepShearTest(0.002);
            break;
          case 2:
            basisBuild.stepShearTest(-0.002);
            break;
          case 3:
            basisBuild.stretch2(-0.002);
            break;
        }
        basisBuild.stepQuasistatic();
        VECTOR displacements = basisBuild.getDisplacement();
        int u_size = displacements.size();
        if( i == 0 && j == 0 ) fprintf(file, "%i %i\n", u_size, 60);

        for(int k = 0; k < u_size; k++)
          fprintf(file, "%lf ", displacements[k]);

        fprintf(file, "\n");
      }
    }
    fclose(file);
    // build the scene
    triangleMesh.buildBlob(sceneNum, argv[1], false, basisCols);
  }
  else
  {
    // build the scene
    triangleMesh.buildBlob(sceneNum, argv[1], true, 0);
  }

  bodyForce[0] = 0;
  bodyForce[1] = -0.3;

  triangleMesh.addWall(WALL(VEC2(1,0), VEC2(-0.98,0)));
  triangleMesh.addWall(WALL(VEC2(-1,0), VEC2(0.98,0)));
  triangleMesh.addWall(WALL(VEC2(0,1), VEC2(0,-0.95)));

}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  cout << " Usage: " << argv[0] << "<vertex filename> <which scene> <extra args>" << endl;
  cout << "\t Valid values: " << endl;
  cout << "\t\t <which test>: SINGLE, LSHEAR, RSHEAR, SQUASH, STRETCH, MOTION" << endl;

  readCommandLine(argc, argv);

  // printf("-----------Finite Difference test on Startup------------\n");
  // HessianDifference();
  // PK1Difference();
  // printf("---------------End Finite Difference test---------------\n");

  // initialize GLUT and GL
  glutInit(&argc, argv);

  // open the GL window
  glutWindow();
  return 1;
}
