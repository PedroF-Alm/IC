#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <chrono>
#include "include/matrix2D.h"

#define N 1000
#define N_A 1000
#define M_A 1000
#define N_B 1000
#define M_B 1000
#define MIN 1
#define MAX 10

#define PARALLEL 1
// #define GENERATE 1
// #define PRINT 1

matrix2D multiply(matrix2D A, matrix2D B);
void generate_mat(matrix2D M, int min, int max);
void print_mat(matrix2D M);

int main()
{
    omp_set_num_threads(4);

    #ifdef N
        matrix2D A = matrix2D(N, N), B = matrix2D(N, N);
    #else
        matrix2D A = matrix2D(N_A, M_A), B = matrix2D(N_B, M_B);
    #endif


    #ifdef GENERATE
        generate_mat(A, MIN, MAX);
        generate_mat(B, MIN, MAX);
    #endif

    auto begin = std::chrono::high_resolution_clock::now();
    matrix2D C = multiply(A, B);
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    printf("t = %fs\n", duration / 1000.0);

    #ifdef PRINT
        print_mat(A);
        print_mat(B);
        print_mat(C);
    #endif

    return 0;
}

matrix2D multiply(matrix2D A, matrix2D B)
{
    matrix2D C = matrix2D(A.get_n(), B.get_m());

    if (A.get_m() != B.get_n())
        return C;

    int i, j, k, aux = 0;
    #ifdef PARALLEL
        #pragma omp parallel for collapse(2) default(none) private(i, j, k) shared(A, B, C) reduction(+:aux)
    #endif
    for (i = 0; i < C.get_n(); i++)
        for (j = 0; j < C.get_m(); j++)
        {
            int aux = 0;
            for (k = 0; k < A.get_m(); k++)                
                aux += A.get(i, k) * B.get(k, j);
            C.set(i, j, aux);
        }

    return C;
}

void generate_mat(matrix2D M, int min, int max)
{
    srand(0);

    int i, j;
    #ifdef PARALLEL
        #pragma omp parallel for collapse(2) private(i, j) shared(M)
    #endif
    for (i = 0; i < M.get_n(); i++)
        for (j = 0; j < M.get_m(); j++)
            M.set(i, j, rand() % (max - min) + min);
}

void print_mat(matrix2D M)
{
    for (int i = 0; i < M.get_n(); i++)
    {
        for (int j = 0; j < M.get_m(); j++)
            printf("%d\t", M.get(i, j));
        printf("\n");
    }
}