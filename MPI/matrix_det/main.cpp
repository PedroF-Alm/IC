#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define N 200
#define MAX 20
#define MIN 0

//#define PRINT

int **matrix = NULL;

int determinant(int**, int, int, int);

int main(int argc, char **argv)
{
    srand(0);

    int min = MIN, max = MAX;

    if (argc > 2)
    {
        min = atoi(argv[1]);
        max = atoi(argv[2]);
    }

    matrix = (int**) calloc(N, sizeof(int*));
    for (int i = 0; i < N; i++)
    {
        matrix[i] = (int*) calloc(N, sizeof(int));
        for (int j = 0; j < N; j++)
            matrix[i][j] = (int)((double)rand() / RAND_MAX * (max - min + 1)) + min;
    }

    // MPI
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_size == 1)
    {
        MPI_Finalize();
        printf("Cannot run on only one process!\n");
        return 0;
    }

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // ____________________________________________

    int det = determinant(matrix, N, world_rank, world_size);
    printf("%d\n", det);
    // ____________________________________________

    MPI_Finalize();

    for (int i = 0; i < N; i++)
        free(matrix[i]);
    free(matrix);

    return 0;
}

// Function to calculate the determinant
// of a matrix
int determinant(int **matrix, int size, int world_rank, int world_size)
{
    int det = 0, aux = 0;
    int sign = 1;

    if (world_rank != 0)
        MPI_Recv(&matrix, size*size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Base Case
    if (size == 1)
    {
        det = matrix[0][0];
    }
    else if (size == 2)
    {
        det = (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
    }

    // Perform the Laplace Expansion
    else
    {
        for (int i = 0; i < size; i++)
        {

            // Stores the cofactor matrix
            int **cofactor = new int *[size - 1];
            for (int j = 0; j < size - 1; j++)
            {
                cofactor[j] = new int[size - 1];
            }
            int sub_i = 0, sub_j = 0;
            for (int j = 1; j < size; j++)
            {
                for (int k = 0; k < size; k++)
                {
                    if (k == i)
                    {
                        continue;
                    }
                    cofactor[sub_i][sub_j] = matrix[j][k];
                    sub_j++;
                }
                sub_i++;
                sub_j = 0;
            }


            if (world_rank == 0) {
                MPI_Send(&cofactor, (size - 1) * (size - 1), MPI_INT, i % world_size + 1, 0, MPI_COMM_WORLD);
            }
            else {
                det += sign * matrix[0][i] * determinant(cofactor, size - 1, world_rank, world_size); 
            }
            sign = -sign;

            for (int j = 0; j < size - 1; j++)
            {
                delete[] cofactor[j];
            }
            delete[] cofactor;
        }

        if (world_rank == 0)
        {
            for (int i = 0; i < size; i++)
            {
                MPI_Recv(&aux, i % world_size - 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                det += aux;
            }
        }
        else
        {
            MPI_Send(&det, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        }
    }

    // Return the final determinant value
    return det;
}