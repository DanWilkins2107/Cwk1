#include "csr_matrix.hpp"
#include "mesh.hpp"
#include "vector.hpp"

csr_matrix CreateMatrixStorageFromMesh(mesh input_mesh)
{
    // Creating csr_matrix
    csr_matrix matrix;

    // Working out how many non-zero entries matrix_a will have.
    int nonzero_entries = 0;
    for (int i = 0; i < input_mesh.no_cells; i++)
    {
        // Add one for self
        nonzero_entries++;
        for (int j = 0; j < 4; j++)
        {
            if (input_mesh.cells[i].neighbours[j] != -1)
            {
                nonzero_entries++;
            }
        }
    }

    // Allocating Memory
    matrix.no_rows = input_mesh.no_cells;
    matrix.matrix_entries = new double[nonzero_entries];
    matrix.column_no = new int[nonzero_entries];
    matrix.row_start = new int[input_mesh.no_cells + 1];

    // Add values to column_no and row_start
    int total_entries = 0;
    matrix.row_start[0] = 0;
    for (int i = 0; i < input_mesh.no_cells; i++)
    {
        matrix.column_no[total_entries] = i;
        total_entries++;
        for (int j = 0; j < 4; j++)
        {
            if (input_mesh.cells[i].neighbours[j] != -1)
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

double* CreateVectorStorageFromMesh(mesh input_mesh)
{
    // Allocating Memory
    double* vec;
    vec = new double[input_mesh.no_cells];

    // Returning the vector
    return vec;
}

csr_matrix PopulateMatrixA(mesh input_mesh, double* vector_f, double (*f)(double, double), double (*g)(double, double), double* (*b)(double, double))
{
    csr_matrix matrix_a = CreateMatrixStorageFromMesh(input_mesh);
    for (int i = 0; i < input_mesh.no_cells; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            double direction_term = ComputeDotProduct(b(input_mesh.cells[i].edge_centres[j][0], input_mesh.cells[i].edge_centres[j][1]), input_mesh.cells[i].edge_normals[j], 2);
            double h_numerator = input_mesh.cells[i].mesh_widths[j % 2];
            double h_denominator = input_mesh.cells[i].mesh_widths[(j + 1) % 2];
            // Note: D in variable names refers to generalised direction
            if (input_mesh.cells[i].neighbours[j] != -1)
            {
                // Case 1: Not a Boundary
                // Note N: j = 0, W: j = 1, S: j = 2: E: j = 3;
                double u_C_coefficient = h_numerator / h_denominator;
                double u_D_coefficient = -1.0 * u_C_coefficient;
                if (direction_term >= 0) {
                    u_C_coefficient += h_numerator * direction_term;
                } else {
                    u_D_coefficient += h_numerator * direction_term;
                }
            }
            else
            {
                // Case 2: A Boundary
            }
        }
    }

    return matrix_a;
}

double* PopulateVectorF(mesh input_mesh, double (*f)(double, double), double (*g)(double, double), double* (*b)(double, double))
{
    double* vector_f = CreateVectorStorageFromMesh(input_mesh);

    for (int i = 0; i < input_mesh.no_cells; i++)
    {
        double h_x = input_mesh.cells[i].mesh_widths[0];
        double h_y = input_mesh.cells[i].mesh_widths[1];
        double xbar_i = input_mesh.cells[i].cell_centre[0];
        double ybar_i = input_mesh.cells[i].cell_centre[1];

        vector_f[i] = h_x * h_y * f(xbar_i, ybar_i);
    }

    return vector_f;
}