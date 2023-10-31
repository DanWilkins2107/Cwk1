//////////////////////////////////////////////////////////////
//Module that implements a sparse matrix diagonal data structure
//////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <math.h>
#include <cassert>
#include "vector.hpp"
#include "dense_matrix.hpp"


//////////////////////////////////////////////////////////////
void AllocateDenseMatrix(dense_matrix& matrix,int no_rows,int no_cols)
//Deallocate storage of sparse matrix and set other entries to 0
//for completeness
{
	matrix.no_rows = no_rows;
	matrix.no_cols = no_cols;

	matrix.matrix_entries = new double*[no_rows];

	for (int k=0;k<matrix.no_rows;k++)
	{
		matrix.matrix_entries[k] = new double[no_cols];

		for (int j=0;j<matrix.no_cols;j++)
		{
			matrix.matrix_entries[k][j] = 0.0;
		}
	}	
}

//////////////////////////////////////////////////////////////
void PrintDenseMatrix(dense_matrix& matrix)
//Prints matrix to the screen in a non-pretty way
{

	for (int k=0;k<matrix.no_rows;k++)
	{
		for (int j=0;j<matrix.no_cols;j++)
		{
			std::cout << matrix.matrix_entries[k][j] << " ";
		}

		std::cout << std::endl;
	}	

	std::cout << std::endl;
}

//////////////////////////////////////////////////////////////
void DeleteDenseMatrix(dense_matrix& matrix)
//Deallocate storage of sparse matrix and set other entries to 0
//for completeness
{
	for (int k=0;k<matrix.no_rows;k++)
	{
		delete[] matrix.matrix_entries[k];
	}
	delete[] matrix.matrix_entries;

	matrix.matrix_entries = NULL;
	matrix.no_rows = 0;
	matrix.no_cols = 0;
}

//////////////////////////////////////////////////////////////
void AddColumnTriangularMatrix(dense_matrix& matrix,double* column)
//Appends an extra column to the right side of an upper triangular
//matrix, thus keepin the matrix upper triangular
{
	int no_rows = matrix.no_rows;
	int no_cols = matrix.no_cols;
	double** matrix_entries = matrix.matrix_entries;

	double** new_matrix_entries = new double*[no_rows+1];

	for (int k=0;k<no_rows;k++)
	{
		double* new_row = AllocateVector(no_cols+1);
		for (int j=0;j<no_cols;j++)
		{
			new_row[j] = matrix_entries[k][j];
		}
		new_row[no_cols] = column[k];
		
		new_matrix_entries[k] = new_row;
	}

	//Last row that is mainly zeros
	double* new_row = AllocateVector(no_cols+1);
	new_row[no_cols] = column[no_rows];

	new_matrix_entries[no_rows] = new_row;

	//Delete the old matrix entries
	if (matrix_entries != NULL)
	{	
		for (int k=0;k<no_rows;k++)
		{
			delete[] matrix_entries[k];
		}

		delete[] matrix_entries;
	}

	//Update with the new matrix entries
	matrix.matrix_entries = new_matrix_entries;
	matrix.no_rows++;
	matrix.no_cols++;
}

//////////////////////////////////////////////////////////////
void BackwardSubstitution(dense_matrix matrix, double* rhs, double* solution)
//Function to solve Ax = b, where A is an upper triangular matrix
//matrix - A
//rhs - B
//x - solution (assumed preallocated)
{
	int no_rows = matrix.no_rows;
	int no_cols = matrix.no_cols;
	double** matrix_entries = matrix.matrix_entries;

	for (int k=no_rows-1;k>=0;k--)
	{
		solution[k] = rhs[k];

		for (int j=no_cols-1;j>k;j--)
		{
			solution[k] -= matrix_entries[k][j]*solution[j];
		}

		solution[k] /= matrix_entries[k][k];
	}
}

//////////////////////////////////////////////////////////////
void MultiplyDenseMatrixVector(dense_matrix matrix, double* vector, double* result)
//Function to multiply a vector from the left by a dense matrix
//Assumes the memory for the resultant vector is preallocated
{
	//Set some new variables to make code easier to read
	int no_rows = matrix.no_rows;
	int no_cols = matrix.no_cols;
	double** matrix_entries = matrix.matrix_entries;

	//Zero result - must be done if memory is reused
	for (int k=0;k<matrix.no_rows;k++)
	{
		result[k] = 0.0;
	}

	for (int k=0;k<no_rows;k++)
	{
		for (int j=0;j<no_cols;j++) //Loop over only correct entries in each diagonal
		{
			result[k] += matrix_entries[k][j]*vector[j];
		}
	}
}

//////////////////////////////////////////////////////////////
void AddRowToMatrix(dense_matrix& matrix,double* new_row, int n)
//Dynamically adds a row (new row) to a dense matrix (matrix)
//n is the length of the new row.
{
	double** matrix_entries = matrix.matrix_entries;
	int no_rows = matrix.no_rows;
	double** updated_matrix = new double*[no_rows+1];

	// If this is an unallocated matrix, we need to 
	// set the number of columns
	if (matrix_entries == NULL)
	{
		matrix.no_cols = n;
	}
	else
	{
		if (n != matrix.no_cols)
		{
			std::cout << "Matrix size mismatch" << std::endl;
			return;
		}
	}

	//Store pointers to the old matrix
	for (int i=0;i<no_rows;i++)
	{
		updated_matrix[i] = matrix_entries[i];
	}

	//Add the new row
	updated_matrix[no_rows] = new_row;

	//Clean up the old storage
	if (matrix_entries != NULL)
	{
		delete[] matrix_entries;
	}

	matrix.matrix_entries = updated_matrix;
	matrix.no_rows++;
}

//////////////////////////////////////////////////////////////
double* PerformGaussianElimination(dense_matrix& matrix, double* RHS)
{
  //Functions to solve the matrix problem Ax = b
	//A is stored in matrix
	//b is stored in RHS

	//We assume that matrix is square and RHS is also of a
	//compatible size
	

	int no_rows = matrix.no_rows;
	double** matrix_vals = matrix.matrix_entries;
  double* solution = AllocateVector(no_rows); //Set up solution vector

  //Perform Elimination in place in matrix - using partial pivoting

  for (int i=0;i<no_rows;i++)
  {
    //Find the maximum element in the column
    
    double max_val = fabs(matrix_vals[i][i]);
    int index = i;
    for (int k=i+1;k<no_rows;k++)
    {
      if (fabs(matrix_vals[k][i]) > max_val)
      {
        max_val = fabs(matrix_vals[k][i]);
        index = k;
      }
    }

    //Swap Rows in matrix
    double* row_ptr_tmp = matrix_vals[index];
    matrix_vals[index] = matrix_vals[i];
    matrix_vals[i] = row_ptr_tmp;

    //Swap rows in the RHS
    double rhs_tmp = RHS[index];
    RHS[index] = RHS[i];
    RHS[i] = rhs_tmp;
    
    //Now do the elimination 
    double multiplier = 1.0/matrix_vals[i][i];

    //Scale the current row
    RHS[i] *= multiplier;
    for (int j=i;j<no_rows;j++)
    {
      matrix_vals[i][j] *= multiplier;
    }

    //Eliminate on remaining rows
    for (int k=i+1;k<no_rows;k++)
    {
      for (int j=i+1;j<no_rows;j++)
      {
        matrix_vals[k][j] -= matrix_vals[k][i]*matrix_vals[i][j]; 
      }
      RHS[k] -= matrix_vals[k][i]*RHS[i];
    }
  }

	//NOw solve the upper trianglar matrix problem
	BackwardSubstitution(matrix, RHS, solution);

  	return solution;

}



