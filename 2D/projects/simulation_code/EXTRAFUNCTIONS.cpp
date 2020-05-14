#include "EXTRAFUNCTIONS.h"
#include <iostream>

//general function to vectorize a matrix. will be used throughout code.
VECTOR vectorize(MATRIX F)
{
  int rows = F.rows();
  int cols = F.cols();
  VECTOR ret( rows * cols);
  int count = 0;

  for(int j = 0; j < cols; j++ )
  {
    for(int i = 0; i < rows; i++ )
    {
      ret[count] = F(i, j);
      count++;
    }
  }
  return ret;
}

///general function to print out a matrix. will be used to debug only
void printMatrix(MATRIX matrix)
{
  int rows = matrix.rows();
  int cols = matrix.cols();
  for(int i = 0; i < rows; i++)
  {
    printf("row %d: ( ", i+1);
    for(int j = 0; j < cols; j++)
    {
      printf("%f, ", matrix(i,j));
    }
    printf(")\n");
  }
}

void printVector(VECTOR vec)
{
  int size = vec.size();
  printf("( ");
  for(int i =0; i < size; i++)
  {
    printf("%f ", vec[i]);
  }
  printf(")\n");
}

//double contraction function for 2 same sized matrices
Real doubleContraction(MATRIX m1, MATRIX m2)
{
  MATRIX c = m1.cwiseProduct(m2);
  return c.sum();
};
