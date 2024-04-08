#include "../include/Mult.h"

int multiply(int *A, int *B, int n) {
    int result = 0;
    for (int i = 0; i < n; i++)
        result += A[i] * B[i];
    return result;
}

void freemat(int **M, int n) {
    for (int i = 0; i < n; i++)
        if (M[i] != NULL)
            free(M[i]), M[i] = NULL;
    free(M);
}

int **allocatemat(int n) {
    int **M = (int**) calloc(n, sizeof(int*));
    for (int i = 0; i < n; i++)
        M[i] = (int*) calloc(n, sizeof(int));
    return M;
}