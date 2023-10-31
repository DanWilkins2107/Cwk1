//Code to show use of Gaussian Elimination

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include "vector.hpp"
#include "dense_matrix.hpp"

int main(int argc, char* argv[])
{
  dense_matrix matrix; //matrix
  double* RHS; // Allocatable RHS vector
  int no_rows = 10;

  //Set up a random matrix and random RHS vector
  AllocateDenseMatrix(matrix,no_rows,no_rows);
  RHS = AllocateVector(no_rows);

  for (int i=0;i<no_rows;i++)
  {
    RHS[i] = ((double) rand() / (RAND_MAX));

    for (int j=0;j<no_rows;j++)
    {
      matrix.matrix_entries[i][j] = ((double) rand() / (RAND_MAX));
    }
  }
  
  //Now apply Gaussian elimination to solve Ax = RHS
  double* solution = PerformGaussianElimination(matrix,RHS);

  PrintVector(solution,no_rows);

  //Now deallocate matrix and vectors
  DeallocateVector(RHS);
  DeallocateVector(solution);
  DeleteDenseMatrix(matrix);

  return 0;
}