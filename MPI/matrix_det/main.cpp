#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define N 3
#define MAX 20
#define MIN 0

// #define PRINT

int **matrix = NULL, **submatrix = NULL;

int determinant(int **, int);

int main(int argc, char **argv)
{
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

    if (world_rank == 0)
    {
        int det = 0;

        matrix = (int **)malloc(sizeof(int *) * N);
        submatrix = (int **)malloc(sizeof(int *) * (N - 1));

        srand(0);

        int min = MIN, max = MAX;

        if (argc > 2)
        {
            min = atoi(argv[1]);
            max = atoi(argv[2]);
        }
        for (int i = 0; i < N; i++)
        {
            matrix[i] = (int *)malloc(sizeof(int) * N);
            if (i < N - 1)
                submatrix[i] = (int *)malloc(sizeof(int) * (N - 1));
            for (int j = 0; j < N; j++)
                matrix[i][j] = (int)((double)rand() / RAND_MAX * (max - min + 1)) + min;
        }

        // expand along first row
        for (int col = 0; col < N; col++)
        {
            // copy into minor matrix, left side then right side
            for (int r = 1; r < N; r++)
            {
                for (int c = 0; c < col; c++)
                {
                    submatrix[r - 1][c] = matrix[r][c];
                }
                for (int c = col + 1; c < N; c++)
                {
                    submatrix[r - 1][c - 1] = matrix[r][c];
                }
            }

            // use "minor" matrix at this point to calculte
            // its determinant
            MPI_Send(&(matrix[0][col]), 1, MPI_INT, 1 + (col % world_size - 1), 0, MPI_COMM_WORLD);
            MPI_Send(submatrix, (N - 1) * (N - 1), MPI_INT, 1 + (col % world_size - 1), 1, MPI_COMM_WORLD);
        }

        for (int worker = 1; worker < world_size; worker++)
        {
            int d;
            MPI_Recv(&d, 1, MPI_INT, worker, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            det += d;
        }

        printf("Determinant: %d\n", det);

        for (int i = 0; i < N; i++)
            free(matrix[i]);
        free(matrix);

        for (int i = 0; i < N - 1; i++)
            free(submatrix[i]);
        free(submatrix);
    }
    else
    {
        submatrix = (int **)malloc(sizeof(int *) * (N - 1));
        for (int i = 0; i < N - 1; i++)
            submatrix[i] = (int *)malloc(sizeof(int) * (N - 1));

        int det = 0;
        for (int col = 0; col < N; col += world_size - 1)
        {
            int cofactor;
            int s = col % 2 == 0 ? -1 : -1;
            MPI_Recv(&cofactor, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(submatrix, (N - 1) * (N - 1), MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            det += s * cofactor * determinant(matrix, N - 1);
        }

        MPI_Send(&det, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);

        for (int i = 0; i < N - 1; i++)
            free(submatrix[i]);
        free(submatrix);
    }

    // ____________________________________________

    MPI_Finalize();

    return 0;
}

// Function to calculate the determinant
// of a matrix
int determinant(int **matrix, int size)
{
    int det = 0;
    int sign = 1;

    if (matrix == NULL)
        return 0;

    for (int i = 0; i < size; i++)
        if (matrix[i] == NULL)
            return 0;

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

            // Update the determinant value
            det += sign * matrix[0][i] * determinant(cofactor, size - 1);
            sign = -sign;
            for (int j = 0; j < size - 1; j++)
            {
                delete[] cofactor[j];
            }
            delete[] cofactor;
        }
    }

    // Return the final determinant value
    return det;
}