#include "csr_matrix.hpp"
#include "mesh.hpp"
#include "finite_volume.hpp"
#include <iostream>

int main () {
    // Ask user for number of cells in each direction.
    int number_of_cells_x, number_of_cells_y;
    std::cout << "How many cells in the x direction?" << std::endl;
    std::cin >> number_of_cells_x;
    std::cout << "How many cells in the y direction?" << std::endl;
    std::cin >> number_of_cells_y;

    // Setting Up Mesh
    mesh rectangular_mesh = ConstructRectangularMesh(
        number_of_cells_x, number_of_cells_y, 0, 0, 0, 0
    );

    //Create matrix storage from mesh
    csr_matrix mesh_matrix = CreateMatrixStorageFromMesh(rectangular_mesh);

    //Print mesh
    PrintCSRMatrix(mesh_matrix);
}