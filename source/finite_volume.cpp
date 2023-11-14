#include "csr_matrix.hpp"
#include "mesh.hpp"

csr_matrix CreateMatrixStorageFromMesh(mesh mesh)
{
    // Creating csr_matrix
    csr_matrix matrix;

    // Working out how many non-zero entries matrix_a will have.
    int nonzero_entries = 0;
    for (int i = 0; i < mesh.no_cells; i++)
    {
        //Add one for self
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
    matrix.no_rows = mesh.no_cells;
    matrix.matrix_entries = new double[nonzero_entries];
    matrix.column_no = new int[nonzero_entries];
    matrix.row_start = new int[mesh.no_cells + 1];

    // Add values to column_no and row_start
    int total_entries = 0;
    matrix.row_start[0] = 0;
    for (int i = 0; i < mesh.no_cells; i++)
    {
        matrix.column_no[total_entries] = i;
        total_entries++;
        for (int j = 0; j < 4; j++)
        {
            if (mesh.cells[i].neighbours[j] != -1)
            {
                matrix.column_no[total_entries] = j;
                total_entries++;
            }
        }
        matrix.row_start[i + 1] = total_entries;
    }
    
    // Returning the matrix
    return matrix;
}

double* CreateVectorStorageFromMesh(mesh mesh)
{
    // Allocating Memory
    double* vec;
    vec = new double[mesh.no_cells];

    // Returning the vector
    return vec;
}
