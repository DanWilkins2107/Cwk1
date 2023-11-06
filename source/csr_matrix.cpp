#include "../include/csr_matrix.hpp"
#include <iostream>

// Function Declarations
csr_matrix SetupMatrixA();
void DeallocateCSRMatrix(csr_matrix matrix);
void MultiplyMatrixVector(csr_matrix& matrix, double* vector, double* productVector);

// Function to set up and return matrix A.
csr_matrix SetupMatrixA()
{
    // Creating matrix_a
    csr_matrix matrix_a;

    // Setting up values
    double matrix_entries[] = {8, 2, 3, 1, 4, 6, 7};
    int column_no[] = {0, 3, 1, 2, 2, 0, 3};
    int row_start[] = {0, 2, 4, 5, 7};

    // Allocating Memory to Arrays
    matrix_a.matrix_entries = new double[7];
    matrix_a.column_no = new int[7];
    matrix_a.row_start = new int[5];

    // Setting Up Values
    for (int i = 0; i < 7; i++)
    {
        matrix_a.matrix_entries[i] = matrix_entries[i];
    }

    for (int i = 0; i < 7; i++)
    {
        matrix_a.column_no[i] = column_no[i];
    }

    for (int i = 0; i < 5; i++)
    {
        matrix_a.row_start[i] = row_start[i];
    }

    matrix_a.no_rows = 4;

    // Returning matrix_a
    return matrix_a;
}

// Function to deallocate memory from any CSR Matrix.
void DeallocateCSRMatrix(csr_matrix matrix)
{
    delete[] matrix.column_no;
    delete[] matrix.matrix_entries;
    delete[] matrix.row_start;
}

// Function to take in a matrix and vector, and change product vector to the product.
void MultiplyMatrixVector(csr_matrix& matrix, double* vector, double* productVector)
{
    // Convert CSR matrix to regular form matrix
    double** regular_matrix;
    regular_matrix = new double*[matrix.no_rows];
    for (int i = 0; i < matrix.no_rows; i++)
    {
        regular_matrix[i] = new double[matrix.no_rows];
        for (int j = matrix.row_start[i]; j < matrix.row_start[i + 1]; j++)
        {
            regular_matrix[i][matrix.column_no[j]] = matrix.matrix_entries[j];
        }
    }

    // Perform Matrix Multiplication
    for (int i = 0; i < matrix.no_rows; i++)
    {
        double count = 0;
        for (int j = 0; j < matrix.no_rows; j++)
        {
            count += regular_matrix[i][j] * vector[j];
        }
        productVector[i] = count;
    }
}
