//////////////////////////////////////////////////////////////
//Module that implements a double precision vector and a 
//number of vector operations
//////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>

//////////////////////////////////////////////////////////////
double* AllocateVector(int n)
//Allocates a double precision vector and sets all entries to zero
{
	double* vector = new double[n];
	for (int i=0;i<n;i++)
	{
		vector[i] = 0.0;
	}

	return vector;
}

//////////////////////////////////////////////////////////////
void ZeroVector(double* vector, int n)
//zeros a vector and sets all entries to zero
{
	for (int i=0;i<n;i++)
	{
		vector[i] = 0.0;
	}
}

//////////////////////////////////////////////////////////////
void ScaleVector(double* vector,double scaling,int n)
//Scales a vector
{
	for (int i=0;i<n;i++)
	{
		vector[i] *= scaling;
	}
}
//////////////////////////////////////////////////////////////
void PrintVector(double* vector,int n)
//Prints a vector to the screen
{
	for (int i=0;i<n;i++)
	{
		std::cout << vector[i] << std::endl;
	}
}

//////////////////////////////////////////////////////////////
void CopyVector(double* vector,double* copied_vector,int n)
//Makes a copy of a vector - assumes that the memory has already
//been allocated
{
	for (int k=0;k<n;k++)
	{
		copied_vector[k] = vector[k];
	}
}
//////////////////////////////////////////////////////////////
void SubtractVectors(double* vec1,double* vec2,int n)
//Overwrites vec1 with vec1-vec2
{
	for (int k=0;k<n;k++)
	{
		vec1[k] -= vec2[k];
	} 
}

//////////////////////////////////////////////////////////////
double ComputeDotProduct(double* vector1,double* vector2,int n)
//Compute the dot product of two vectors
{
	double dot_prod = 0.0; 
	for (int k=0;k<n;k++)
	{
		dot_prod += vector1[k]*vector2[k];
	}

	return dot_prod;
}

//////////////////////////////////////////////////////////////
void CombineVectors(double* vector1,double* vector2,double scale,int n)
//Carries out the operation x -> x+scale*y, where x = vector1, y = vector2
{
	for (int k=0;k<n;k++)
	{
		vector1[k] += scale*vector2[k];
	}
}

//////////////////////////////////////////////////////////////
void MultiplyVectorsElementwise(double* vec1,double* vec2,int n)
//Overwrites vec1 with vec1*vec2 elementwise
{
	for (int k=0;k<n;k++)
	{
		vec1[k] *= vec2[k];
	}
}

//////////////////////////////////////////////////////////////
double NormVector(double* vector,int n)
//Compute the l2-norm of a vector
{
	double norm = 0.0;
	for (int k=0;k<n;k++)
	{
		norm += pow(vector[k],2);
	}

	return sqrt(norm);
}

//////////////////////////////////////////////////////////////
double FindMaximum(double* vector,int n)
//Function to find maximum of a vector
{
	double max_val = vector[0];
	for (int k=1;k<n;k++)
	{
		if (vector[k] > max_val)
		{   
			max_val = vector[k];
		}
	}

	return max_val;
}

//////////////////////////////////////////////////////////////
void DeallocateVector(double* vector)
//Deletes storage for a vector
{
	delete[] vector;
}

//////////////////////////////////////////////////////////////
double* ReadVectorFromFile(std::string fileName)
//Read a vector from file
{  
	std::ifstream read_file;
	read_file.open(fileName); // Open file
	assert(read_file.is_open()); 

	std::string line;

	//Read no of rows
	std::getline(read_file, line); //Ignore String
	std::getline(read_file, line);
	int n = atoi(line.c_str());

	double* vector = AllocateVector(n);

	std::getline(read_file, line); //Ignore String
	for (int k=0;k<n;k++)
	{
		std::getline(read_file, line);
		vector[k] = atof(line.c_str());
	}

	// Finally close file
	read_file.close();

	return vector;
}