#ifndef csr_matrix_header
#define csr_matrix_header
#include <string>

struct csr_matrix
{
    // CSR format matrix
	double* matrix_entries;
    int* column_no;
    int* row_start;

    // Explicitly defined number of rows 
    int no_rows;
};

void MultiplyMatrixVector(csr_matrix& matrix, double* vector, double* productVector);
void DeallocateCSRMatrix(csr_matrix matrix);
csr_matrix ReadMatrix(std::string matrix_filename);
double* ReadVector(std::string vector_filename);
csr_matrix SetupMatrixA();

#endif