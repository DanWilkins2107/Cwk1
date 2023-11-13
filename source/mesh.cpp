//////////////////////////////////////////////////////////////
//Module that implements a double precision matrix and a 
//number of matrix operations
//////////////////////////////////////////////////////////////

#include "dense_matrix.hpp"
#include "mesh.hpp"

//////////////////////////////////////////////////////////////
mesh ConstructRectangularMesh(int Nx, int Ny, double a, double b, 
	double c, double d)
// Function to construct a rectangular mesh on the region [axb]x[cxd]
//Nx - Number of cells in x-direction. Must be even and >= 2
//Ny - Number of cells in y-direction. Must be even and >= 2
{
	mesh mesh_data;
	
	int number_nodes = (Nx+1)*(Ny+1);
	int number_cells = Nx*Ny;

	mesh_data.no_nodes = number_nodes;
	mesh_data.no_cells = number_cells;

	//First set up nodes
	
	AllocateDenseMatrix(mesh_data.coordinates,number_nodes,2);

	double hx = (b-a)/((double) Nx);
	double hy = (d-c)/((double) Ny);

	int count = 0;

	//Nodes ordered left to right, then down to up
	for (int j=0;j<Ny+1;j++)
	{
		double y_coord = c+j*hy;

		for (int i=0;i<Nx+1;i++)
		{
			double x_coord = a+i*hx;

			mesh_data.coordinates.matrix_entries[count][0] = x_coord;
			mesh_data.coordinates.matrix_entries[count][1] = y_coord;

			count++;
		}
	}

	//Now set up cells
	//Again ordered left to right, then down to up
	mesh_data.cells = new cell_information[number_cells];
	count = 0;

	int SW = 0;
	int SE = 1;
	int NW = Nx+1;
	int NE = Nx+2;

	int E, W, S, N;

	for (int j=0;j<Ny;j++)
	{
		for (int i=0;i<Nx;i++)
		{
			//Set up nodes
			mesh_data.cells[count].nodes = new int[4];
			mesh_data.cells[count].nodes[0] = NE++;
			mesh_data.cells[count].nodes[1] = NW++;
			mesh_data.cells[count].nodes[2] = SW++;
			mesh_data.cells[count].nodes[3] = SE++;

			//Set up neighbours
			mesh_data.cells[count].neighbours = new int[4];
			mesh_data.cells[count].edge_normals = new double*[4];
			mesh_data.cells[count].edge_centres = new double*[4];

			if (i==0) //Left Boundary
			{
				W = -1;
				E = count+1;
			}
			else if (i==Nx-1) //Right Boundary
			{
				W = count-1;
				E = -1;
			}
			else //Interior
			{
				W = count-1;
				E = count+1;
			}

			if (j==0) //Lower Boundary
			{
				S = -1;
				N = count+Nx;
			}
			else if (j==Ny-1) //Upper Boundary
			{
				S = count-Nx;
				N = -1;
			}
			else
			{
				S = count-Nx;
				N = count+Nx;
			}

			mesh_data.cells[count].neighbours[0] = N;
			mesh_data.cells[count].neighbours[1] = W;
			mesh_data.cells[count].neighbours[2] = S;
			mesh_data.cells[count].neighbours[3] = E;

			// Mesh widths
			mesh_data.cells[count].mesh_widths = new double[2];
			mesh_data.cells[count].mesh_widths[0] = hx;
			mesh_data.cells[count].mesh_widths[1] = hy;

			//Cell centre
			mesh_data.cells[count].cell_centre = new double[2];
			mesh_data.cells[count].cell_centre[0] = hx/((double)2)+i*hx;
			mesh_data.cells[count].cell_centre[1] = hy/((double)2)+j*hy;

			//Edge centres and edge normals

			SetEdgeInformation(mesh_data.cells[count],hx,hy);
		
			count++;
		}
		NE++;
		NW++;
		SW++;
		SE++;
	}

	return mesh_data;
}



mesh ConstructLShapedMesh(int Nx, int Ny, double a, double b, double c, double d)
// Function to construct an L-shaped mesh on the region [axb]x[cxd], with the 
// top right quadrant removed
//Nx - Number of cells in x-direction. Must be even and >= 4
//Ny - Number of cells in y-direction. Must be even and >= 4
{
	mesh mesh_data;
	
	int number_nodes = (Nx+1)*(Ny+1)-(Nx*Ny)/4;
	int number_cells = Nx*Ny-(Nx*Ny)/4;

	mesh_data.no_nodes = number_nodes;
	mesh_data.no_cells = number_cells;

	//First set up nodes
	
	AllocateDenseMatrix(mesh_data.coordinates,number_nodes,2);

	double hx = (b-a)/((double) Nx);
	double hy = (d-c)/((double) Ny);

	int count = 0;

	//Nodes ordered left to right, then down to up

	// Lower rectangular portion
	for (int j=0;j<Ny/2+1;j++)
	{
	double y_coord = c+j*hy;

	for (int i=0;i<Nx+1;i++)
	{
		double x_coord = a+i*hx;

		mesh_data.coordinates.matrix_entries[count][0] = x_coord;
		mesh_data.coordinates.matrix_entries[count][1] = y_coord;

		count++;
	}
	}

	//Top left portion
	for (int j=Ny/2+1;j<Ny+1;j++)
	{
	double y_coord = c+j*hy;

	for (int i=0;i<Nx/2+1;i++)
	{
		double x_coord = a+i*hx;

		mesh_data.coordinates.matrix_entries[count][0] = x_coord;
		mesh_data.coordinates.matrix_entries[count][1] = y_coord;

		count++;
	}
	}

	//Now set up cells
	//Again ordered left to right, then down to up
	mesh_data.cells = new cell_information[number_cells];
	count = 0;

	// Lower Portion, excluding top row

	int SW = 0;
	int SE = 1;
	int NW = Nx+1;
	int NE = Nx+2;

	int E, W, S, N;

	for (int j=0;j<Ny/2-1;j++)
	{
	for (int i=0;i<Nx;i++)
	{
		//Set up nodes
		mesh_data.cells[count].nodes = new int[4];
		mesh_data.cells[count].nodes[0] = NE++;
		mesh_data.cells[count].nodes[1] = NW++;
		mesh_data.cells[count].nodes[2] = SW++;
		mesh_data.cells[count].nodes[3] = SE++;

		//Set up neighbours
		mesh_data.cells[count].neighbours = new int[4];
		mesh_data.cells[count].edge_normals = new double*[4];
		mesh_data.cells[count].edge_centres = new double*[4];

		if (i==0) //Left Boundary
		{
			W = -1;
			E = count+1;
		}
		else if (i==Nx-1) //Right Boundary
		{
			W = count-1;
			E = -1;
		}
		else //Interior
		{
			W = count-1;
			E = count+1;
		}

		if (j==0) //Lower Boundary
		{
			S = -1;
			N = count+Nx;
		}
		else if (j==Ny-1) //Upper Boundary
		{
			S = count-Nx;
			N = -1;
		}
		else
		{
			S = count-Nx;
			N = count+Nx;
		}


		mesh_data.cells[count].neighbours[0] = N;
		mesh_data.cells[count].neighbours[1] = W;
		mesh_data.cells[count].neighbours[2] = S;
		mesh_data.cells[count].neighbours[3] = E;

		// Mesh widths
		mesh_data.cells[count].mesh_widths = new double[2];
		mesh_data.cells[count].mesh_widths[0] = hx;
		mesh_data.cells[count].mesh_widths[1] = hy;

		//Cell centre
		mesh_data.cells[count].cell_centre = new double[2];
		mesh_data.cells[count].cell_centre[0] = hx/((double)2)+i*hx;
		mesh_data.cells[count].cell_centre[1] = hy/((double)2)+j*hy;

		//Edge centres and edge normals

		SetEdgeInformation(mesh_data.cells[count],hx,hy);

		count++;
	}
	NE++;
	NW++;
	SW++;
	SE++;
	}

	// Connecting row 1

	for (int i=0;i<Nx;i++)
	{
		//Set up nodes
		mesh_data.cells[count].nodes = new int[4];
		mesh_data.cells[count].nodes[0] = NE++;
		mesh_data.cells[count].nodes[1] = NW++;
		mesh_data.cells[count].nodes[2] = SW++;
		mesh_data.cells[count].nodes[3] = SE++;

		//Set up neighbours
		mesh_data.cells[count].neighbours = new int[4];
		mesh_data.cells[count].edge_normals = new double*[4];
		mesh_data.cells[count].edge_centres = new double*[4];

		if (i==0) //Left Boundary
		{
			W = -1;
			E = count+1;
		}
		else if (i==Nx-1) //Right Boundary
		{
			W = count-1;
			E = -1;
		}
		else //Interior
		{
			W = count-1;
			E = count+1;
		}

		if (i >= Nx/2) //Right boundary
		{
			S = std::max(count-Nx,-1);
			N = -1;
		}
		else
		{
			S = std::max(count-Nx,-1);
			N = count+Nx;
		}


		mesh_data.cells[count].neighbours[0] = N;
		mesh_data.cells[count].neighbours[1] = W;
		mesh_data.cells[count].neighbours[2] = S;
		mesh_data.cells[count].neighbours[3] = E;

		// Mesh widths
		mesh_data.cells[count].mesh_widths = new double[2];
		mesh_data.cells[count].mesh_widths[0] = hx;
		mesh_data.cells[count].mesh_widths[1] = hy;

		//Cell centre
		mesh_data.cells[count].cell_centre = new double[2];
		mesh_data.cells[count].cell_centre[0] = hx/((double)2)+i*hx;
		mesh_data.cells[count].cell_centre[1] = hy/((double)2)+(Ny/2-1)*hy;

		//Edge centres and edge normals
		SetEdgeInformation(mesh_data.cells[count],hx,hy);

		count++;
	}
	NE++;
	NW++;
	SW++;
	SE++;

	//Top left portion 

	// Connecting row 2

	for (int i=0;i<Nx/2;i++)
	{
		//Set up nodes
		mesh_data.cells[count].nodes = new int[4];
		mesh_data.cells[count].nodes[0] = NE++;
		mesh_data.cells[count].nodes[1] = NW++;
		mesh_data.cells[count].nodes[2] = SW++;
		mesh_data.cells[count].nodes[3] = SE++;

		//Set up neighbours
		mesh_data.cells[count].neighbours = new int[4];
		mesh_data.cells[count].edge_normals = new double*[4];
		mesh_data.cells[count].edge_centres = new double*[4];

		if (i==0) //Left Boundary
		{
			W = -1;
			E = count+1;
		}
		else if (i==Nx/2-1) //Right Boundary
		{
			W = count-1;
			E = -1;
		}
		else //Interior
		{
			W = count-1;
			E = count+1;
		}

		S = count-Nx;
		N = count+Nx/2;
		
		mesh_data.cells[count].neighbours[0] = N;
		mesh_data.cells[count].neighbours[1] = W;
		mesh_data.cells[count].neighbours[2] = S;
		mesh_data.cells[count].neighbours[3] = E;

		// Mesh widths
		mesh_data.cells[count].mesh_widths = new double[2];
		mesh_data.cells[count].mesh_widths[0] = hx;
		mesh_data.cells[count].mesh_widths[1] = hy;

		//Cell centre
		mesh_data.cells[count].cell_centre = new double[2];
		mesh_data.cells[count].cell_centre[0] = hx/((double)2)+i*hx;
		mesh_data.cells[count].cell_centre[1] = hy/((double)2)+(Ny/2)*hy;

		//Edge centres and edge normals
		SetEdgeInformation(mesh_data.cells[count],hx,hy);

		count++;
	}
	NE++;
	NW++;
	SW += Nx/2+1;
	SE += Nx/2+1;

	//Other cells

	for (int j=Ny/2+1;j<Ny;j++)
	{
	for (int i=0;i<Nx/2;i++)
	{
		//Set up nodes
		mesh_data.cells[count].nodes = new int[4];
		mesh_data.cells[count].nodes[0] = NE++;
		mesh_data.cells[count].nodes[1] = NW++;
		mesh_data.cells[count].nodes[2] = SW++;
		mesh_data.cells[count].nodes[3] = SE++;

		//Set up neighbours
		mesh_data.cells[count].neighbours = new int[4];
		mesh_data.cells[count].edge_normals = new double*[4];
		mesh_data.cells[count].edge_centres = new double*[4];

		if (i==0) //Left Boundary
		{
			W = -1;
			E = count+1;
		}
		else if (i==Nx/2-1) //Right Boundary
		{
			W = count-1;
			E = -1;
		}
		else //Interior
		{
			W = count-1;
			E = count+1;
		}
		
		if (j==Ny-1) //Upper Boundary
		{
			S = count-Nx/2;
			N = -1;
		}
		else
		{
			S = count-Nx/2;
			N = count+Nx/2;
		}


		mesh_data.cells[count].neighbours[0] = N;
		mesh_data.cells[count].neighbours[1] = W;
		mesh_data.cells[count].neighbours[2] = S;
		mesh_data.cells[count].neighbours[3] = E;

		// Mesh widths
		mesh_data.cells[count].mesh_widths = new double[2];
		mesh_data.cells[count].mesh_widths[0] = hx;
		mesh_data.cells[count].mesh_widths[1] = hy;

		//Cell centre
		mesh_data.cells[count].cell_centre = new double[2];
		mesh_data.cells[count].cell_centre[0] = hx/((double)2)+i*hx;
		mesh_data.cells[count].cell_centre[1] = hy/((double)2)+j*hy;

		//Edge centres and edge normals
		SetEdgeInformation(mesh_data.cells[count],hx,hy);

		count++;
	}
	NE++;
	NW++;
	SW++;
	SE++;
	}

	return mesh_data;
}

//////////////////////////////////////////////////////////////
void SetEdgeInformation(cell_information& cell,double hx, double hy)
//Helper function to set up the edge information given a cell
{

	// Neighbour 1
		cell.edge_centres[0] = new double[2];
		cell.edge_centres[0][0] = 
		cell.cell_centre[0];
		cell.edge_centres[0][1] = 
		cell.cell_centre[1]+hy/2.0;

		cell.edge_normals[0] = new double[2];
		cell.edge_normals[0][0] = 0.0;
		cell.edge_normals[0][1] = 1.0;

		// Neighbour 2
		cell.edge_centres[1] = new double[2];
		cell.edge_centres[1][0] = 
		cell.cell_centre[0]-hx/2.0;
		cell.edge_centres[1][1] = 
		cell.cell_centre[1];

		cell.edge_normals[1] = new double[2];
		cell.edge_normals[1][0] = -1.0;
		cell.edge_normals[1][1] = 0.0;

		// Neighbour 3
		cell.edge_centres[2] = new double[2];
		cell.edge_centres[2][0] = 
		cell.cell_centre[0];
		cell.edge_centres[2][1] = 
		cell.cell_centre[1]-hy/2.0;

		cell.edge_normals[2] = new double[2];
		cell.edge_normals[2][0] = 0.0;
		cell.edge_normals[2][1] = -1.0;

		// Neighbour 3
		cell.edge_centres[3] = new double[2];
		cell.edge_centres[3][0] = 
		cell.cell_centre[0]+hx/2.0;
		cell.edge_centres[3][1] = 
		cell.cell_centre[1];

		cell.edge_normals[3] = new double[2];
		cell.edge_normals[3][0] = 1.0;
		cell.edge_normals[3][1] = 0.0;
}

//////////////////////////////////////////////////////////////
void DeallocateMesh(mesh meshData)
//Deletes/resets all data in the mesh structure
{
	for (int k=0;k<meshData.no_cells;k++)
	{
	delete[] meshData.cells[k].nodes;
	delete[] meshData.cells[k].neighbours;
	delete[] meshData.cells[k].mesh_widths;

	meshData.cells[k].nodes = NULL;
	meshData.cells[k].neighbours = NULL;
	meshData.cells[k].mesh_widths = NULL;

	for (int j=0;j<4;j++)
	{
		delete[] meshData.cells[k].edge_centres[j];
		delete[] meshData.cells[k].edge_normals[j];
	}

	delete[] meshData.cells[k].edge_centres;
	delete[] meshData.cells[k].edge_normals;

	meshData.cells[k].edge_centres = NULL;
	meshData.cells[k].edge_normals = NULL;
	}

	DeleteDenseMatrix(meshData.coordinates);
	delete[] meshData.cells;
	meshData.cells = NULL;

	meshData.no_cells = 0;
	meshData.no_nodes = 0;
}

//////////////////////////////////////////////////////////////
void PrintMesh(mesh meshData)
//Print some of the mesh data structure to the screen.
{
	std::cout << "No nodes = " << meshData.no_nodes << std::endl;
	std::cout << "No cells = " << meshData.no_cells << std::endl;

	std::cout << "\nNode coordinates" << std::endl;
	PrintDenseMatrix(meshData.coordinates);
	
	std::cout << std::endl;

	//Now print a selection of the cell information
	for (int k=0;k<meshData.no_cells;k++)
	{
	std::cout << "Cell = " << k << std::endl;
	std::cout << "Nodes = " << meshData.cells[k].nodes[0] << " "
		<< " " << meshData.cells[k].nodes[1] << " " << meshData.cells[k].nodes[2]
		<< " " << meshData.cells[k].nodes[3] << std::endl << std::endl;
	std::cout << "Neighbours = " << meshData.cells[k].neighbours[0] << " "
		<< meshData.cells[k].neighbours[1] << " " << meshData.cells[k].neighbours[2]
		<< " " << meshData.cells[k].neighbours[3] << std::endl << std::endl;
	}

}

void OutputSolution(double* solution, mesh mesh_data, std::string filename)
{
		//Function to output a solution vector as matrix for reading into
		//Matlab or Python

		std::ofstream outputFile;
		outputFile.open(filename);

    int no_cells = mesh_data.no_cells;

		for (int k=0;k < no_cells;k++)
		{
				for (int l=0;l < 4; l++)
				{
						outputFile << mesh_data.coordinates.matrix_entries[mesh_data.cells[k].nodes[l]][0] << " ";
            outputFile << mesh_data.coordinates.matrix_entries[mesh_data.cells[k].nodes[l]][1] << " ";
        }

        outputFile << solution[k];

				outputFile << std::endl;
		}

		outputFile.close();

}