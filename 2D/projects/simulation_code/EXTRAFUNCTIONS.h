#ifndef EXTRAFUNCTIONS_H
#define EXTRAFUNCTIONS_H

#include "SETTINGS.h"

// vectorize function
VECTOR vectorize(const MATRIX F);

//print out general matrix
void printMatrix(MATRIX matrix);

void printVector(VECTOR vec);
//double contraction function for 2 same sized matrices
Real doubleContraction(MATRIX m1, MATRIX m2);


#endif
