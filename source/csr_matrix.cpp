#include "csr_matrix.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

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

csr_matrix ReadMatrix(std::string matrix_filename)
{
    // Open the File
    std::ifstream read_file(matrix_filename);
    assert(read_file.is_open());

    // Find the number of rows and the number of nonzeros
    std::string current_string;
    int no_rows;
    int no_nonzeros;
    read_file >> current_string >> current_string >> no_rows;
    read_file >> current_string >> current_string >> no_nonzeros;

    // Create CSR Matrix
    csr_matrix matrix;
    matrix.no_rows = no_rows;
    matrix.row_start = new int[no_rows + 1];
    matrix.column_no = new int[no_nonzeros];
    matrix.matrix_entries = new double[no_nonzeros];

    // Add values to row_start
    read_file >> current_string;
    for (int i = 0; i < no_rows + 1; i++)
    {
        read_file >> matrix.row_start[i];
    }

    // Add values to column_no
    read_file >> current_string;
    for (int i = 0; i < no_nonzeros; i++)
    {
        read_file >> matrix.column_no[i];
    }

    // Add values to matrix_entries
    read_file >> current_string;
    for (int i = 0; i < no_nonzeros; i++)
    {
        read_file >> matrix.matrix_entries[i];
    }

    // Close File
    read_file.close();
    return matrix;
}

double* ReadVector(std::string vector_filename)
{
    // Open the file
    std::ifstream read_file(vector_filename);
    assert(read_file.is_open());

    // Find vector length
    std::string current_string;
    int vector_length;
    read_file >> current_string >> current_string >> vector_length;

    // Allocate memory to new vector
    double* vector;
    vector = new double[vector_length];

    // Fill vector
    read_file >> current_string >> current_string;
    for (int i = 0; i < vector_length; i++)
    {
        read_file >> vector[i];
    }

    return vector;
}

void PrintCSRMatrix(csr_matrix matrix) {
    // Print number of rows
    std::cout << "Number Rows" << std::endl;
    std::cout << matrix.no_rows << std::endl;

    // Print number of nonzero entries in matrix
    std::cout << "Number Nonzeros" << std::endl;

    int no_nonzeros = matrix.row_start[matrix.no_rows];
    std::cout << no_nonzeros << std::endl;
    
    // Print the row start values
    std::cout << "row_start" << std::endl;
    for (int i = 0; i < matrix.no_rows + 1; i++) {
        std::cout << matrix.row_start[i] << std::endl;
    }
    
    // Print the column number values
    std::cout << "column_no" << std::endl;
    for (int i = 0; i < no_nonzeros; i++) {
        std::cout << matrix.column_no[i] << std::endl;
    }

    // Print the nonzero entries
    std::cout << "matrix_entries" << std::endl;
    for (int i = 0; i < no_nonzeros; i++) {
        std::cout << matrix.matrix_entries[i] << std::endl;
    }
}
