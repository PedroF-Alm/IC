#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <chrono>
#include <random>
#include <thread>
#include "include/matrix2D.h"

#define N 3
#define MIN 1
#define MAX 2

#define PARALLEL 4
#define GENERATE 1
#define PRINT 1

int determinant(matrix2D *M);
void generate_mat(matrix2D M, int min, int max);
void print_mat(matrix2D M);

int main()
{
    #ifdef PARALLEL 
        omp_set_num_threads(PARALLEL);
    #endif

    matrix2D M = matrix2D(N, N);

    #ifdef GENERATE
        generate_mat(M, MIN, MAX);        
    #endif

    auto begin = std::chrono::high_resolution_clock::now();
    int det = determinant(&M);
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    printf("t = %fs\n", duration / 1000.0);
    printf("det(M) = %d\n", det);

    #ifdef PRINT
        print_mat(M);
    #endif

    return 0;
}

int determinant(matrix2D *M)
{
    int det = 0;
    int sign = 1;
    int size = M->get_n();
 
    // Base Case
    if (size == 1) {
        det = M->get(0, 0);
    }
    else if (size == 2) {
        det = (M->get(0, 0) * M->get(1, 1))
            - (M->get(0, 1) * M->get(1, 0));
    }
 
    // Perform the Laplace Expansion
    else {
        int i, j, k;
        #ifdef PARALLEL
            #pragma omp parallel for private(i, j, k) shared(size, sign, M) reduction(+: det)
        #endif
        for (int i = 0; i < size; i++) {
            // Stores the cofactor matrix
            matrix2D *cofactor = new matrix2D(size - 1, size - 1);
            int sub_i = 0, sub_j = 0;
            for (j = 1; j < size; j++) {

                for (k = 0; k < size; k++) {
                    if (k == i) {
                        continue;
                    }
                    cofactor->set(sub_i, sub_j, M->get(j, k));
                    sub_j++;
                }
                sub_i++;
                sub_j = 0;
            }
 
            // Update the determinant value
            det += sign * M->get(0, i) * determinant(cofactor);
            sign = -sign;

            delete cofactor;
        }
    }
 
    // Return the final determinant value
    return det;
}

int intRand(const int & min, const int & max) {
    static thread_local std::mt19937 generator;
    std::uniform_int_distribution<int> distribution(min,max);
    return distribution(generator);
}

void generate_mat(matrix2D M, int min, int max)
{
    int i, j;
    #ifdef PARALLEL
        #pragma omp parallel for collapse(2) private(i, j) shared(M)
    #endif
    for (i = 0; i < M.get_n(); i++)
        for (j = 0; j < M.get_m(); j++)
            M.set(i, j, rand() % (max - min + 1) + min);
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