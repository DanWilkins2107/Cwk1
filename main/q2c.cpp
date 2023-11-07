#include "../source/csr_matrix.cpp"

int main() {
    csr_matrix matrix;
    matrix = ReadMatrix("../matrix1.dat");
    double* vector;
    vector = ReadVector("../vector1.dat");
}