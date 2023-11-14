#ifndef finite_volume_header
#define finite_volume_header

#include "csr_matrix.hpp"
#include "mesh.hpp"

csr_matrix CreateMatrixStorageFromMesh(mesh mesh);
double* CreateVectorStorageFromMesh(mesh mesh);
csr_matrix PopulateMatrixA(mesh input_mesh, double (*f)(double, double), double (*g)(double, double), double* (*b)(double, double))
double* PopulateVectorF(mesh input_mesh, double (*f)(double, double), double (*g)(double, double), double* (*b)(double, double))

#endif