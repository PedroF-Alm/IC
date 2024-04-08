#include <stdio.h>
#include <time.h>
#include "./include/Mult.h"

#define N 200
#define MAX 20
#define MIN 0

//#define PRINT

int A[N][N] = {0}, B[N][N] = {0}, C[N][N] = {0};

void master(int);
void worker(int, int);

int main(int argc, char **argv)
{

    srand(0);

    int min = MIN, max = MAX;

    if (argc > 2)
    {
        min = atoi(argv[1]);
        max = atoi(argv[2]);
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i][j] = (int)((double)rand() / RAND_MAX * (max - min + 1)) + min;
            B[i][j] = (int)((double)rand() / RAND_MAX * (max - min + 1)) + min;
        }
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

    if (world_rank == 0)
    {        
        master(world_size);
    }
    else
    {
        worker(world_rank, world_size);
    }

    // ____________________________________________

    MPI_Finalize();

    return 0;
}

void master(int world_size)
{
    printf("Running on %d processes...\n", world_size);
    
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            MPI_Send(&i, 1, MPI_INT, 1 + i % (world_size - 1), 1, MPI_COMM_WORLD);
            MPI_Send(&j, 1, MPI_INT, 1 + i % (world_size - 1), 2, MPI_COMM_WORLD);            
        }
    }    

    for (int i = 0; i < N; i++)
    {
        MPI_Recv(C[i], N, MPI_INT, 1 + i % (world_size - 1), 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    #ifdef PRINT
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                printf("A[%d][%d] = %d\tB[%d][%d] = %d\tC[%d][%d] = %d\n", i, j, A[i][j], j, i, B[i][j], i, j, C[i][j]);
    #endif
}

void worker(int world_rank, int world_size)
{
    int x, y;

    for (int i = world_rank - 1; i < N; i += world_size - 1)
    {        
        for (int j = 0; j < N; j++)
        {
            MPI_Recv(&x, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&y, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            C[i][j] = multiply(A[x], B[y], N);            
        }
    }

    for (int i = world_rank - 1; i < N; i += world_size - 1)
        MPI_Send(C[i], N, MPI_INT, 0, 3, MPI_COMM_WORLD); 
}