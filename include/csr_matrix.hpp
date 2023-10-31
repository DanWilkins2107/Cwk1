#ifndef csr_matrix_setup
#define csr_matrix_setup

struct csr_matrix
{
    // CSR format matrix
	double* matrix_entries;
    double* column_no;
    double* row_start;
    // Explicitly defined number of rows 
    int number_of_rows;
};

#endif