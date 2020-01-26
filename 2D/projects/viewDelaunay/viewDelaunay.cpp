#include "SETTINGS.h"
#include <cmath>
#include <iostream>

#include <vector>

#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
//#include <GLUT/glut.h>
#include <GL/freeglut.h>
#endif

using namespace std;

// the resolution of the OpenGL window -- independent of the field resolution
int xScreenRes = 850;
int yScreenRes = 850;

// Text for the title bar of the window
string windowLabel("Delaunay Viewer");

// mouse tracking variables
int xMouse         = -1;
int yMouse         = -1;
int mouseButton    = -1;
int mouseState     = -1;
int mouseModifiers = -1;

// current grid cell the mouse is pointing at
int xField = -1;
int yField = -1;

// animate the current runEverytime()?
bool animate = false;

// the current viewer eye position
//VEC3 eyeCenter(0.5, 0.5, 1);
//VEC3 eyeCenter(0.2, 0.8, 1);
VEC3 eyeCenter(0, 0, 1);

// current zoom level into the field
float zoom = 1.25;

// the actual triangle data
vector<VEC2> nodes;
vector<int> indices;
vector<VEC3I> triangles;

VEC2 mins, maxs;

///////////////////////////////////////////////////////////////////////
// Print a string to the GL window
///////////////////////////////////////////////////////////////////////
void printGlString(string output)
{
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
  for (unsigned int x = 0; x < output.size(); x++)
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, output[x]);
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

  glColor4f(1,0,0,1);
  glPointSize(10.0); 
  glBegin(GL_POINTS);
    for (unsigned int x = 0; x < nodes.size(); x++)
      glVertex2f(nodes[x][0], nodes[x][1]);
  glEnd();

  glColor4f(0,0,1,1);
  glBegin(GL_TRIANGLES);
    for (unsigned int x = 0; x < triangles.size(); x++)
    {
      int node0 = triangles[x][0];
      int node1 = triangles[x][1];
      int node2 = triangles[x][2];
      glVertex2f(nodes[node0][0], nodes[node0][1]);
      glVertex2f(nodes[node1][0], nodes[node1][1]);
      glVertex2f(nodes[node2][0], nodes[node2][1]);
    }
  glEnd();

  glutSwapBuffers();
}

///////////////////////////////////////////////////////////////////////
// Map the arrow keys to something here
///////////////////////////////////////////////////////////////////////
void glutSpecial(int key, int x, int y)
{
  switch (key)
  {
    case GLUT_KEY_LEFT:
        break;
    case GLUT_KEY_RIGHT:
        break;
    case GLUT_KEY_UP:
        break;
    case GLUT_KEY_DOWN:
        break;
    default:
        break;
  }
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
  
  float xDiff = x - xMouse;
  float yDiff = y - yMouse;
  float speed = 0.001;
  
  if (mouseButton == GLUT_LEFT_BUTTON) 
  {
    eyeCenter[0] -= xDiff * speed;
    eyeCenter[1] += yDiff * speed;
  }
  if (mouseButton == GLUT_RIGHT_BUTTON)
    zoom -= yDiff * speed;
  
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
  }
  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glutWindow()
{
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE| GLUT_RGBA);
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
  glutSpecialFunc(&glutSpecial);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);
  glutPassiveMotionFunc(&glutPassiveMouseMotion);
 
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  // enter the infinite GL loop
  glutMainLoop();
  
  // Control flow will never reach here
  return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////
// Read in a *.node file
///////////////////////////////////////////////////////////////////////
void readNodes(const string& filename, vector<VEC2>& nodes, vector<int>& indices)
{
  // read the nodes file
  FILE* file = NULL;
  file = fopen(filename.c_str(), "r");

  if (file == NULL)
  {
    cout << " File " << filename.c_str() << " does not exist! " << endl;
    exit(0);
  }

  // read in the total number of nodes and other attributes
  int totalNodes = -1;
  int dimension = -1;
  int totalAttributes = -1;
  int totalBoundaryMarkers = -1;
  fscanf(file, "%i %i %i %i", &totalNodes, &dimension, &totalAttributes, &totalBoundaryMarkers);

  cout << " Total nodes: " << totalNodes << endl;
  cout << " Dimension: " << dimension << endl;
  cout << " Attributes: " << totalAttributes << endl;
  cout << " Boundary markers: " << totalBoundaryMarkers << endl;

  assert(dimension == 2);

  for (int x = 0; x < totalNodes; x++)
  {
    // get the vertex position
    int index = -1;
    double position[2];
    fscanf(file, "%i %lf %lf", &index, &(position[0]), &(position[1]));
    cout << " index: " << index << "\t node: " << position[0] << " " << position[1] << endl;

    // store it as a node
    VEC2 node(position[0], position[1]);
    nodes.push_back(node);
    indices.push_back(index);

    // strip off the attributes
    double throwAway;
    for (int y = 0; y < totalAttributes; y++)
      fscanf(file, "%lf", &throwAway);

    // strip off the boundary markers
    for (int y = 0; y < totalBoundaryMarkers; y++)
      fscanf(file, "%lf", &throwAway);
  }
  fclose(file);
}

///////////////////////////////////////////////////////////////////////
// Read in a *.ele file
///////////////////////////////////////////////////////////////////////
void readElements(const string& filename, const int offset, vector<VEC3I>& triangles)
{
  FILE* file = NULL;
  file = fopen(filename.c_str(), "r");

  if (file == NULL)
  {
    cout << " File " << filename.c_str() << " does not exist! " << endl;
    exit(0);
  }

  int totalTriangles = -1;
  int totalNodesPerTriangle = -1;
  int totalAttributes = -1;
  fscanf(file, "%i %i %i", &totalTriangles, &totalNodesPerTriangle, &totalAttributes);

  cout << " Total triangles: " << totalTriangles << endl;
  cout << " Total nodes in each triangle: " << totalNodesPerTriangle << endl;
  cout << " Total attributes per triangle: " << totalAttributes << endl;
  assert(totalNodesPerTriangle == 3);

  for (int x = 0; x < totalTriangles; x++)
  {
    int triangleIndex = -1;
    int nodeIndices[3];

    fscanf(file, "%i %i %i %i", &triangleIndex, &nodeIndices[0], &nodeIndices[1], &nodeIndices[2]);

    VEC3I triangle(nodeIndices[0], nodeIndices[1], nodeIndices[2]);
    triangle -= VEC3I(offset, offset, offset);
    triangles.push_back(triangle);

    cout << " Triangle " << triangleIndex << ": " << triangle[0] << " " << triangle[1] << " " << triangle[2] << endl;
  }
}

///////////////////////////////////////////////////////////////////////
// Load up a 2D triangle mesh
///////////////////////////////////////////////////////////////////////
void loadTriangles2D(const string& prefix)
{
  // read the nodes file
  string nodeFile = prefix + string(".node");
  readNodes(nodeFile, nodes, indices);

  // did the file indices start at 1 or 0? If 1, we need to subtract off one
  // from all the triangles' node indices.
  int offset = (indices[0] == 0) ? 0 : 1;

  // read in the elements file
  string elementFile = prefix + string(".ele");
  readElements(elementFile, offset, triangles);
}

///////////////////////////////////////////////////////////////////////
// Get the bounding box of the nodes
///////////////////////////////////////////////////////////////////////
void getBoundingBox(VEC2 mins, VEC2 maxs)
{
  mins = nodes[0];
  maxs = nodes[0];

  for (unsigned int x = 1; x < nodes.size(); x++)
  {
    if (nodes[x][0] < mins[0])
      mins[0] = nodes[x][0];
    if (nodes[x][1] < mins[1])
      mins[1] = nodes[x][1];
    if (nodes[x][0] > maxs[0])
      maxs[0] = nodes[x][0];
    if (nodes[x][1] > maxs[1])
      maxs[1] = nodes[x][1];
  }

  cout << " Bounding box mins: " << mins[0] << " " << mins[1] << endl;
  cout << " Bounding box maxs: " << maxs[0] << " " << maxs[1] << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  if (argc < 2)
  {
    cout << " USAGE: " << argv[0] << " <Triangle output prefix> " << endl;
    return 0;
  }

  loadTriangles2D(argv[1]);
  getBoundingBox(mins, maxs);

  eyeCenter[0] = (mins[0] + maxs[0]) * 0.5;
  eyeCenter[1] = (mins[1] + maxs[1]) * 0.5;

  // initialize GLUT and GL
  glutInit(&argc, argv);
  
  // open the GL window
  glutWindow();
  return 1;
}
