#ifndef linear_algebra_header
#define linear_algebra_header

#include "vector.hpp"
#include "dense_matrix.hpp"

//////////////////////////////////////////////////////////////
//Prototypes
//////////////////////////////////////////////////////////////
void PerformArnoldiIteration(csr_matrix& matrix, dense_matrix& krylov_matrix, int k, double* hessenberg);
void ComputeGivensRotations(dense_matrix& givensRotations, double rho, double sigma);
void ApplyGivensRotations(dense_matrix& givensRotations,double* vector);
double* PerformGMRESRestarted(csr_matrix& matrix, double* rhsVector, double* x0,
	int max_iterations, double tol, int restart=100);

#endif
