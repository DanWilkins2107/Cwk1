#include "csr_matrix.hpp"

// Function Declarations
csr_matrix SetupMatrixA();
void DeallocateCSR_Matrix(csr_matrix matrix);
void MultiplyMatrixVector(csr_matrix& matrix, double* vector, double* productVector);

// Function to set up and return matrix A.
csr_matrix SetupMatrixA()
{
    csr_matrix matrix_a;
    matrix_a.column_no = {8, 2, 3, 1, 4, 6, 7};
    matrix_a.matrix_entries = {0, 3, 1, 2, 2, 0, 3};
    matrix_a.row_start = {0, 2, 4, 5, 7} matrix_a.number_of_rows = 4;

    return matrix_a;
}

// Function to deallocate memory from any CSR Matrix.
void DeallocateCSR_Matrix(csr_matrix matrix)
{
    delete[] matrix.column_no;
    delete[] matrix.matrix_entries;
    delete[] matrix.row_start;
}

// Function to take in a matrix and vector, and change product vector to the product.
void MultiplyMatrixVector(csr_matrix& matrix, double* vector, double* productVector)
{
    
}
