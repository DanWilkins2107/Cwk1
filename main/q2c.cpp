#include "../source/linear_algebra.cpp"
#include <iostream>

int main()
{
    // Read in A and x
    csr_matrix matrix_a;
    matrix_a = ReadMatrix("../matrix2.dat");
    double* vector_x;
    vector_x = ReadVector("../vector2.dat");

    // Find b
    double* vector_b;
    vector_b = new double[64];
    MultiplyMatrixVector(matrix_a, vector_x, vector_b);

    // Create initial guess vector x_0
    double* x_0;
    x_0 = new double[64];
    for (int i = 0; i < 64; i++) {
        x_0[i] = 0.5;
    }

    // Perform GMRES algorithm
    double* approximation;
    approximation = PerformGMRESRestarted(matrix_a, vector_b, x_0);

    // Find error and print
    SubtractVectors(vector_x, approximation, 64);
    double error;
    error = NormVector(vector_x, 64);
    std::cout << "The error is " << error;

    // Deallocate Memory
    DeallocateCSRMatrix(matrix_a);
    DeallocateVector(vector_x);
    DeallocateVector(vector_b);
    DeallocateVector(x_0);
    DeallocateVector(approximation);
}