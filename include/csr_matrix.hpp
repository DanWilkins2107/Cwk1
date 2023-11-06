#ifndef csr_matrix_setup
#define csr_matrix_setup

struct csr_matrix
{
    // CSR format matrix
	double* matrix_entries;
    int* column_no;
    int* row_start;

    // Explicitly defined number of rows 
    int no_rows;
};

#endif