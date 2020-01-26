#LDFLAGS_COMMON = -framework Accelerate -framework GLUT -framework OpenGL -L/usr/local/lib/x86_64 -L/opt/local/lib -lstdc++ -lz -fopenmp -msse2
LDFLAGS_COMMON = -framework Accelerate -framework GLUT -framework OpenGL -L/usr/local/lib/x86_64 -L/opt/local/lib -lstdc++ -lz -msse2
#MY_INCLUDES = -I../../src/linearalgebra -I./ -I/opt/local/include -I../../ -I../../src/ -I../../src/solvers/ -I../../src/geometry/ -I../../src/materials -I../../src/util/ -I../../src/solvers/ -I/opt/X11/include
MY_INCLUDES = -I../../src/linearalgebra -I./ -I/opt/local/include -I../../ -I../../src/ -I../../src/solvers/ -I../../src/geometry/ -I../../src/materials -I../../src/util/ -I../../src/solvers/ -I/opt/X11/include -I../../src/spectra-0.5.0/include/

# fixes the inline exceeded warnings
#CFLAGS_COMMON = -c -Wall -I../../src/linearalgebra -I./ -I/opt/local/include -I../../ -O3 -fopenmp -msse2
CFLAGS_COMMON = -c -Wall -I../../src/linearalgebra -I./ -I/opt/local/include -I../../ -O3 -msse2
#COMPILER = clang
COMPILER = g++
