#ifndef dense_matrix_header
#define dense_matrix_header



#include <iostream>
#include <fstream>
#include <math.h>
#include "vector.hpp"

//////////////////////////////////////////////////////////////
struct dense_matrix
//Structure to store a sparse matrix in diagonal form.
//The matrix is assumed to be square
{
	int no_rows = 0; //Number of rows in the matrix.
	int no_cols = 0; //Number of columns in the matrix
	double** matrix_entries = NULL; //All the matrix entries
};

//////////////////////////////////////////////////////////////
//Function Prototypes
//////////////////////////////////////////////////////////////
void AllocateDenseMatrix(dense_matrix& matrix,int no_rows,int no_cols);
void DeleteDenseMatrix(dense_matrix& matrix);
void PrintDenseMatrix(dense_matrix& matrix);
void MultiplyDenseMatrixVector(dense_matrix matrix, double* vector, double* result);
void AddColumnTriangularMatrix(dense_matrix& matrix,double* column);
void BackwardSubstitution(dense_matrix matrix, double* rhs, double* solution);
void AddRowToMatrix(dense_matrix& matrix,double* new_row, int n);
double* PerformGaussianElimination(dense_matrix& matrix, double* RHS);

#endif