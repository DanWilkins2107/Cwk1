#include "csr_matrix.hpp"
#include "mesh.hpp"

csr_matrix CreateMatrixStorageFromMesh(mesh mesh);
double* CreateVectorStorageFromMesh(mesh mesh);

csr_matrix CreateMatrixStorageFromMesh(mesh mesh)
{
    // Creating csr_matrix
    csr_matrix matrix_a;

    // Working out how many non-zero entries matrix_a will have.
    int nonzero_entries = 0;
    for (int i = 0; i < mesh.no_cells; i++)
    {
        nonzero_entries++;
        for (int j = 0; j < 4; j++)
        {
            if (mesh.cells[i].neighbours[j] != -1)
            {
                nonzero_entries++;
            }
        }
    }

    // Allocating Memory
    matrix_a.no_rows = mesh.no_cells;
    matrix_a.matrix_entries = new double[nonzero_entries];
    matrix_a.column_no = new int[nonzero_entries];
    matrix_a.row_start = new int[mesh.no_cells + 1];

    // Returning the matrix
    return matrix_a;
}

double* CreateVectorStorageFromMesh(mesh mesh)
{
    // Allocating Memory
    double* vec;
    vec = new double[mesh.no_cells];

    // Returning the vector
    return vec;
}
