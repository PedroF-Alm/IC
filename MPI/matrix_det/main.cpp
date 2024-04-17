#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define N 15
#define MAX 3
#define MIN 1

// #define PRINT

int *matrix = NULL, *submatrix = NULL;

int determinant(int *, int);
int get(int *, int, int, int);
void set(int *, int, int, int, int);
void set(int *, int, int);

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

    int min = MIN, max = MAX, n = N;

    if (argc > 2)
    {
        n = atoi(argv[1]);
        min = atoi(argv[2]);
        max = atoi(argv[3]);
    }
    else if (argc > 1)
    {
        n = atoi(argv[1]);
    }

    // ____________________________________________

    if (world_rank == 0)
    {
        int det = 0;

        matrix = (int *)malloc(sizeof(int *) * n * n);
        submatrix = (int *)malloc(sizeof(int *) * (n - 1) * (n - 1) + 1);

        srand(0);

        printf("Matrix\n");
        for (int i = 0; i < n * n; i++)
        {
            set(matrix, i, (int)((double)rand() / RAND_MAX * (max - min + 1)) + min);
            printf("%d\t", matrix[i]);
            if ((i + 1) % n == 0)
                printf("\n");
        }

        if (n == 1)
        {
            det = get(matrix, n, 0, 0);
        }
        else if (n == 2)
        {
            det = (get(matrix, n, 0, 0) * get(matrix, n, 1, 1)) - (get(matrix, n, 0, 1) * get(matrix, n, 1, 0));
        }
        else
        {
            // expand along first row
            for (int col = 0; col < n; col++)
            {
                // copy into minor matrix, left side then right side
                for (int r = 1; r < n; r++)
                {
                    for (int c = 0; c < col; c++)
                    {
                        set(submatrix, n - 1, r - 1, c, get(matrix, n, r, c));
                    }
                    for (int c = col + 1; c < n; c++)
                    {
                        set(submatrix, n - 1, r - 1, c - 1, get(matrix, n, r, c));
                    }
                }

                // use "minor" matrix at this point to calculte
                // its determinant
                submatrix[(n - 1) * (n - 1)] = matrix[col] * (col % 2 == 0 ? 1 : -1);                
                MPI_Send(&(submatrix[0]), (n - 1) * (n - 1) + 1, MPI_INT, 1 + col % (world_size - 1), 1, MPI_COMM_WORLD);
            }            

            for (int worker = 1; worker < world_size; worker++)
            {
                int d;
                MPI_Recv(&d, 1, MPI_INT, worker, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                det += d;
            }
        }

        printf("Determinant: %d\n", det);

        free(matrix);
        free(submatrix);
    }
    else if (n > 2)
    {
        submatrix = (int *)malloc(sizeof(int *) * (n - 1) * (n - 1) + 1);

        int det = 0;
        
        for (int col = 0; col < n; col += world_size - 1)
        {
            MPI_Recv(&(submatrix[0]), (n - 1) * (n - 1) + 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            det += submatrix[(n - 1) * (n - 1)] * determinant(submatrix, n - 1);            
        }
        
        MPI_Send(&det, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        free(submatrix);
    }

    // ____________________________________________

    MPI_Finalize();

    return 0;
}

// Function to calculate the determinant
// of a matrix
int determinant(int *matrix, int size)
{
    int det = 0;
    int sign = 1;

    if (matrix == NULL)
        return 0;

    // Base Case
    if (size == 1)
    {
        det = get(matrix, size, 0, 0);
    }
    else if (size == 2)
    {
        det = (get(matrix, size, 0, 0) * get(matrix, size, 1, 1)) - (get(matrix, size, 0, 1) * get(matrix, size, 1, 0));
    }

    // Perform the Laplace Expansion
    else
    {
        for (int i = 0; i < size; i++)
        {

            // Stores the cofactor matrix
            int *cofactor = (int *)malloc(sizeof(int *) * (size - 1) * (size - 1));
            int sub_i = 0, sub_j = 0;
            for (int j = 1; j < size; j++)
            {
                for (int k = 0; k < size; k++)
                {
                    if (k == i)
                    {
                        continue;
                    }
                    set(cofactor, size - 1, sub_i, sub_j, get(matrix, size, j, k));
                    sub_j++;
                }
                sub_i++;
                sub_j = 0;
            }

            // Update the determinant value
            det += sign * get(matrix, size, 0, i) * determinant(cofactor, size - 1);
            sign = -sign;

            free(cofactor);
        }
    }

    // Return the final determinant value
    return det;
}

int get(int *matrix, int n, int i, int j)
{
    return matrix[i * n + j];
}

void set(int *matrix, int n, int i, int j, int value)
{
    matrix[i * n + j] = value;
}

void set(int *matrix, int i, int value)
{
    matrix[i] = value;
}