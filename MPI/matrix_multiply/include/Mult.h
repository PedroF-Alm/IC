#include <mpi.h>
#include <stdlib.h>

int multiply(int *A, int *B, int n);
void freemat(int **M, int n);
int **allocatemat(int n);