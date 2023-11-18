#include <cassert>

//////////////////////////////////////////////////////////////
// Module that implements the GMRES algorithm
//////////////////////////////////////////////////////////////

#include "csr_matrix.hpp"
#include "dense_matrix.hpp"
#include "linear_algebra.hpp"
#include "vector.hpp"

//////////////////////////////////////////////////////////////

void PerformArnoldiIteration(csr_matrix& matrix, dense_matrix& krylov_matrix, int k, double* hessenberg)
{
    // Check k > 0
    assert(k > 0);

    // Allocate Memory and find the last row of the krylov_matrix
    double* q_k;
    q_k = new double[krylov_matrix.no_cols];
    for (int i = 0; i < krylov_matrix.no_cols; i++)
    {
        q_k[i] = krylov_matrix.matrix_entries[krylov_matrix.no_rows - 1][i];
    } 

    // Allocate Memory and find q_kplus1
    double* q_kplus1;
    q_kplus1 = new double[krylov_matrix.no_cols];
    MultiplyMatrixVector(matrix, q_k, q_kplus1);

    // Fill hessenberg array and edit q_kplus1
    for (int i = 0; i < k; i++)
    {
        hessenberg[i] = ComputeDotProduct(q_kplus1, q_k, krylov_matrix.no_cols);
        CombineVectors(q_kplus1, krylov_matrix.matrix_entries[i], hessenberg[i] * -1.0, krylov_matrix.no_cols);
    }

	// Find h_k+1 and perform scaling
    hessenberg[k] = NormVector(q_kplus1, krylov_matrix.no_cols);
    ScaleVector(q_kplus1, 1.0 / hessenberg[k], krylov_matrix.no_cols);

	// Add row to krylov_matrix
	AddRowToMatrix(krylov_matrix, q_kplus1, krylov_matrix.no_cols);
}

//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
double* PerformGMRESRestarted(csr_matrix& matrix, double* rhsVector, double* x0,
                              int max_iterations, double tol, int restart)
// Carry out the restarted GMRES iteration to solve Ax = b
// matrix - A
// rhsVector - b
// x0 - initial guess vector
// max_iterations - maximum number of iterations to perform
// tol - Stop when ||Ax-b|| < tol
// restart - numner of iteration to perfrom before restart
{
    int n = matrix.no_rows;
    double* residual_vector = AllocateVector(n);
    double* solution = AllocateVector(n);

    CopyVector(x0, solution, n);

    MultiplyMatrixVector(matrix, solution, residual_vector);

    SubtractVectors(residual_vector, rhsVector, n);

    ScaleVector(residual_vector, -1.0, n);

    double error = NormVector(residual_vector, n);

    // std::cout << "||r_0|| = " << error << std::endl;

    dense_matrix krylov_matrix;
    dense_matrix givens_matrix;
    dense_matrix Rmatrix;

    int iteration_count = 0;
    while (iteration_count <= max_iterations)
    {
        double beta = error;
        AllocateDenseMatrix(krylov_matrix, 1, n);
        CopyVector(residual_vector, krylov_matrix.matrix_entries[0], n);
        ScaleVector(krylov_matrix.matrix_entries[0], 1.0 / beta, n);

        for (int k = 1; k <= restart; k++)
        {

            if (error < tol)
            {
                break;
            }
            double* hessenberg = AllocateVector(k + 1);

            PerformArnoldiIteration(matrix, krylov_matrix,
                                    k, hessenberg);

            ApplyGivensRotations(givens_matrix, hessenberg);

            double rho = hessenberg[k - 1];
            double sigma = hessenberg[k];

            ComputeGivensRotations(givens_matrix, rho, sigma);

            hessenberg[k - 1] = sqrt(pow(rho, 2) + pow(sigma, 2));
            hessenberg[k] = 0.0;

            AddColumnTriangularMatrix(Rmatrix, hessenberg);

            // Set up RHS vector to solve
            double* lsq_rhs_vector = AllocateVector(k + 1);

            lsq_rhs_vector[0] = beta;

            ApplyGivensRotations(givens_matrix, lsq_rhs_vector);

            // Set up y-solution and solve
            double* y_solution = AllocateVector(k);

            BackwardSubstitution(Rmatrix, lsq_rhs_vector, y_solution);

            CopyVector(x0, solution, n);
            // Update the solution
            for (int j = 0; j < k; j++)
            {
                CombineVectors(solution, krylov_matrix.matrix_entries[j], y_solution[j], n);
            }

            // Compute the new residual and error
            MultiplyMatrixVector(matrix, solution, residual_vector);

            SubtractVectors(residual_vector, rhsVector, n);

            error = NormVector(residual_vector, n);

            DeallocateVector(hessenberg);

            DeallocateVector(lsq_rhs_vector);

            DeallocateVector(y_solution);

            iteration_count++;
            std::cout << "||r_" << iteration_count << "|| = " << error << std::endl;
        }

        DeleteDenseMatrix(krylov_matrix);
        DeleteDenseMatrix(givens_matrix);
        DeleteDenseMatrix(Rmatrix);

        ScaleVector(residual_vector, -1.0, n);
        // Resest initial guess
        CopyVector(solution, x0, n);

        if (error < tol)
        {
            break;
        }
    }

    DeallocateVector(residual_vector);
    return solution;
}

//////////////////////////////////////////////////////////////
void ComputeGivensRotations(dense_matrix& givensRotations, double rho,
                            double sigma)
// Function to compute the required Given's rotation
{
    double* new_row = new double[2];

    new_row[0] = rho / sqrt(pow(rho, 2) + pow(sigma, 2));
    new_row[1] = sigma / sqrt(pow(rho, 2) + pow(sigma, 2));

    AddRowToMatrix(givensRotations, new_row, 2);
}

//////////////////////////////////////////////////////////////
void ApplyGivensRotations(dense_matrix& givensRotations, double* vector)
// Function to apply the Given's rotations in order to a vector
{
    int number_rotations = givensRotations.no_rows;
    double** matrix_entries = givensRotations.matrix_entries;
    double temp_val1, temp_val2;

    for (int k = 0; k < number_rotations; k++)
    {
        double ck = matrix_entries[k][0];
        double sk = matrix_entries[k][1];

        temp_val1 = ck * vector[k] + sk * vector[k + 1];
        temp_val2 = ck * vector[k + 1] - sk * vector[k];

        vector[k] = temp_val1;
        vector[k + 1] = temp_val2;
    }
}