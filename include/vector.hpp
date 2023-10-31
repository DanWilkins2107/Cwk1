#ifndef vector_header
#define vector_header

#include <iostream>
#include <math.h>
#include <fstream>

//////////////////////////////////////////////////////////////
//Function Prototypes
//////////////////////////////////////////////////////////////

double* AllocateVector(int n);
double* PrintVector(double* vector,int n);
double* ScaleVector(double* vector,double scaling,int n);
void CopyVector(double* vector,double* copied_vector,int n);
void SubtractVectors(double* vec1,double* vec2,int n);
void CombineVectors(double* vector1,double* vector2,double scale,int n);
void MultiplyVectorsElementwise(double* vec1,double* vec2,int n);
double NormVector(double* vector,int n);
double FindMaximum(double* vector,int n);
void DeallocateVector(double* vector);
double* ReadVectorFromFile(std::string fileName);
double ComputeDotProduct(double* vector1,double* vector2,int n);
void ZeroVector(double* vector, int n);

#endif