#include "csr_matrix.hpp"
#include "mesh.hpp"
#include "vector.hpp"
#include "finite_volume.hpp"

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
        int matrix_entries_added = 0;
        for (int j = 0; j < 4; j++)
        {
            double x_n = input_mesh.cells[i].edge_centres[j][0];
            double y_n = input_mesh.cells[i].edge_centres[j][1];

            double direction_term = ComputeDotProduct(b(x_n, y_n), input_mesh.cells[i].edge_normals[j], 2);

            double h_numerator = input_mesh.cells[i].mesh_widths[j % 2];
            double h_denominator = input_mesh.cells[i].mesh_widths[(j + 1) % 2];

            double u_C_coefficient = 0;
            double u_D_coefficient = 0;
            double u_O_coefficient = 0;
            double constant_value = 0;

            double diffusion_sign = (j > 1) ? -1 : 1;
            // D in variable names refers to direction of edge
            // Note: O in variable names refers to opposite direction
            // Note N: j = 0, W: j = 1, S: j = 2: E: j = 3;
            if (input_mesh.cells[i].neighbours[j] != -1)
            {
                // Case 1: Not a Boundary

                // Diffusive terms
                u_C_coefficient += diffusion_sign * h_numerator / h_denominator;
                u_D_coefficient -= u_C_coefficient;

                // Advective terms
                if (direction_term >= 0)
                {
                    u_C_coefficient += h_numerator * direction_term;
                }
                else
                {
                    u_D_coefficient += h_numerator * direction_term;
                }

                // Add values from current edge
                int number_to_add = 0;
                for (int k = 0; k < j; k++)
                {
                    if (input_mesh.cells[i].neighbours[k] != -1)
                    {
                        number_to_add++;
                    }
                }
                matrix_a.matrix_entries[matrix_a.row_start[i] + number_to_add + 1] += u_D_coefficient;
            }
            else
            {
                // Case 2: A Boundary

                // Diffusive terms
                u_C_coefficient -= diffusion_sign * 3.0 * h_numerator / h_denominator;
                u_O_coefficient += diffusion_sign * h_numerator / (3.0 * h_denominator);
                constant_value += diffusion_sign * h_numerator * 8.0 * g(x_n, y_n) / (h_denominator * 3.0);

                // Advective terms
                if (direction_term >= 0)
                {
                    u_C_coefficient += h_numerator * direction_term;
                }
                else
                {
                    constant_value += h_numerator * direction_term * g(x_n, y_n);
                }

                // Add values from opposite edge
                int number_to_add = 0;
                for (int k = 0; k < ((j + 2) % 4); k++)
                {
                    if (input_mesh.cells[i].neighbours[k] != -1)
                    {
                        number_to_add++;
                    }
                }
                matrix_a.matrix_entries[matrix_a.row_start[i] + number_to_add + 1] += u_O_coefficient;

                // Add values to vector_f
                vector_f[i] -= constant_value;
            }

            // Add values to diagonal of matrix a
            matrix_a.matrix_entries[matrix_a.row_start[i]] += u_C_coefficient;
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