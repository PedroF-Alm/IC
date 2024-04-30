#include <stdio.h>
#include <omp.h>

#define N 10000

// #define PRINT

int A[N][N] = {0}, B[N][N] = {0}, C[N][N] = {0};

int main()
{
    omp_set_num_threads(1);

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = 1, B[i][j] = 1, C[i][j] = 1;

    int i, j, k, aux = 0;
    
    #pragma omp parallel for private(i, j, k)     
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
        {
            int aux = 0;
            for (k = 0; k < N; k++)                
                aux += A[i][k] * B[k][j];
            C[i][j] = aux;
        }

    #ifdef PRINT
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                printf("A[%d][%d] = %d\tB[%d][%d] = %d\tC[%d][%d] = %d\n", i, j, A[i][j], j, i, B[i][j], i, j, C[i][j]);
    #endif

    return 0;
}