#include "csr_matrix.hpp"
#include "finite_volume.hpp"
#include "mesh.hpp"

#include <cassert>
#include <iostream>

int main()
{
    // Ask user for number of cells in each direction.
    int number_of_cells_x, number_of_cells_y;
    std::cout << "How many cells in the x direction?" << std::endl;
    std::cin >> number_of_cells_x;
    std::cout << "How many cells in the y direction?" << std::endl;
    std::cin >> number_of_cells_y;

    // Assert no cells even and >= 2
    assert(number_of_cells_x % 2 == 0);
    assert(number_of_cells_y % 2 == 0);
    assert(number_of_cells_x >= 2);
    assert(number_of_cells_y >= 2);

    // Setting up Mmsh
    mesh rectangular_mesh = ConstructRectangularMesh(
        number_of_cells_x, number_of_cells_y, 0, 0, 0, 0);

    // Create matrix storage from mesh
    csr_matrix mesh_matrix = CreateMatrixStorageFromMesh(rectangular_mesh);

    // Print mesh
    PrintCSRMatrix(mesh_matrix);

    // Deallocate memory
    DeallocateCSRMatrix(mesh_matrix);
    DeallocateMesh(rectangular_mesh);
}