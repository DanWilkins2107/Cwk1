#include "csr_matrix.hpp"
#include <iostream>

int main()
{
    // Run function to set up matrix A
    csr_matrix matrix_a = SetupMatrixA();

    // Setup vector x as given in question
    double x[] = {4, -1, 3, 6};
    double* vec_x;
    vec_x = new double[4];
    for (int i = 0; i < 4; i++)
    {
        vec_x[i] = x[i];
    }

    // Setup product vector
    double* product_Ax;
    product_Ax = new double[4];

    // Run multiplication function
    MultiplyMatrixVector(matrix_a, vec_x, product_Ax);

    // Print out resulting product vector
    std::cout << "Ax = " << std::endl;
    for (int i = 0; i < 4; i++)
    {
        std::cout << product_Ax[i] << std::endl;
    }

    // Deallocate memory
    DeallocateCSRMatrix(matrix_a);
    delete[] vec_x;
    delete[] product_Ax;
}