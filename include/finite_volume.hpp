#ifndef finite_volume_header
#define finite_volume_header

#include "csr_matrix.hpp"
#include "mesh.hpp"

csr_matrix CreateMatrixStorageFromMesh(mesh mesh);
double* CreateVectorStorageFromMesh(mesh mesh);

#endif