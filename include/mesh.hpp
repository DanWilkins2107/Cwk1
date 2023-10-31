#ifndef mesh_header
#define mesh_header

#include <iostream>
#include "vector.hpp"
#include "dense_matrix.hpp"

struct cell_information
{
	int* nodes = NULL; //index of nodes associated with the cell,
										//stored in order NE, NW, SW, SE
	int* neighbours = NULL; //Will be a vector storing element neighbours,
													//stored in order N, W, S , W

	double* cell_centre; //Centre of the cell - vector length 2
	double** edge_centres; //Centre of Edges - 4 x 2 matrix
	double** edge_normals; //Edge_normals - 4 x 2 matrix
	double* mesh_widths; //mesh widths (hx, hy) in each coordinate direction
											 //vector length 2
};

struct mesh
{
	int no_nodes = 0; //Number of nodes in the mesh
	int no_cells = 0; //Number of cells in the mesh
	dense_matrix coordinates; //Coordinates of all the nodes
	cell_information* cells; //Information about each cell
};


mesh ConstructRectangularMesh(int Nx, int Ny, double a, double b, double c, double d);
mesh ConstructLShapedMesh(int Nx, int Ny, double a, double b, double c, double d);
void PrintMesh(mesh meshData);
void DeallocateMesh(mesh meshData);
void SetEdgeInformation(cell_information& cell,double hx, double hy);

#endif