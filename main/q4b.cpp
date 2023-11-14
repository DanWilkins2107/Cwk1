#include "csr_matrix.hpp"
#include "finite_volume.hpp"
#include "mesh.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

// All functions needed for Q4b
double return_1(double x, double y)
{
    return 1;
}

double return_0(double x, double y)
{
    return 0;
}

double f2(double x, double y)
{
    return 8.0 * pow(M_PI, 2) * cos(2.0 * M_PI * x) * cos(2.0 * M_PI * y);
}

double g2(double x, double y)
{
    return cos(2.0 * M_PI * x) * cos(2.0 * M_PI * y);
}

double g4(double x, double y)
{
    if (x == 0 && 0.25 < y < 0.75)
    {
        return 1;
    }
    return 0;
}

double* return_0_array(double x, double y)
{
    double* array;
    array = new double[2];
    array[0] = 0;
    array[1] = 0;
    return array;
}

double* return_10_array(double x, double y)
{
    double* array;
    array = new double[2];
    array[0] = 10;
    array[1] = 10;
    return array;
}

double* b4(double x, double y)
{
    double* array;
    array = new double[2];
    array[0] = (50.0 * y) / sqrt(pow(x, 2) + pow(y, 2));
    array[1] = (-50.0 * x) / sqrt(pow(x, 2) + pow(y, 2));
    return array;
}

int main()
{
    // Ask user to choose problem number
    int problem_no;
    std::cout << "What problem number would you like to choose?" << std::endl;
    std::cin >> problem_no;

    // Check the problem number is valid
    assert(problem_no > 0 && problem_no < 5);

    // Ask user to choose number of cells in each direction
    int number_of_cells_x, number_of_cells_y;
    std::cout << "How many cells in the x direction?" << std::endl;
    std::cin >> number_of_cells_x;
    std::cout << "How many cells in the y direction?" << std::endl;
    std::cin >> number_of_cells_y;

    // Initialize problem
    mesh chosen_mesh;

    double (*f)(double, double);
    double (*g)(double, double);
    double* (*b)(double, double);

    switch (problem_no)
    {
    case 1:
        chosen_mesh = ConstructRectangularMesh(number_of_cells_x, number_of_cells_y, 0, 1, 0, 1);
        f = return_1;
        g = return_0;
        b = return_0_array;
        break;
    case 2:
        chosen_mesh = ConstructLShapedMesh(number_of_cells_x, number_of_cells_y, 0, 1, 0, 1);
        f = f2;
        g = g2;
        b = return_0_array;
        break;
    case 3:
        chosen_mesh = ConstructRectangularMesh(number_of_cells_x, number_of_cells_y, -1, 1, -1, 1);
        f = return_1;
        g = return_0;
        b = return_10_array;
        break;
    case 4:
        chosen_mesh = ConstructLShapedMesh(number_of_cells_x, number_of_cells_y, 0, 1, 0, 1);
        f = return_0;
        g = g4;
        b = b4;
        break;
    default:
        break;
    }

    // Find matrices A and F
    csr_matrix matrix_a = PopulateMatrixA(chosen_mesh, f, g, b);
    double* vector_f = PopulateVectorF(chosen_mesh, f, g, b);
}