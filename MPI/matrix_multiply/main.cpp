#include <stdio.h>
#include <mpi.h>

#define N 10000

// #define PRINT

int A[N][N] = {0}, B[N][N] = {0}, C[N][N] = {0};

int main()
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = 1, B[i][j] = 1, C[i][j] = 1;

    // MPI
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // if (world_size == 1)
    // {
    //     MPI_Finalize();
    //     printf("Cannot run on only one process!\n");
    //     return 0;
    // }

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // ____________________________________________

    if (world_rank == 0)
    {        
        printf("Running on %d processes...\n", world_size);

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                // int aux[2] = {i, j};
                // MPI_Send(&aux, 2, MPI_INT, 1 + i % (world_size - 1), 1, MPI_COMM_WORLD);        
                C[i][j] = 0;
                for (int k = 0; k < N; k++)
                    C[i][j] += A[i][k] * B[j][k];
            }
        }    

        // for (int i = 0; i < N; i++)
        // {
        //     MPI_Recv(C[i], N, MPI_INT, 1 + i % (world_size - 1), 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // }

        #ifdef PRINT
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    printf("A[%d][%d] = %d\tB[%d][%d] = %d\tC[%d][%d] = %d\n", i, j, A[i][j], j, i, B[i][j], i, j, C[i][j]);
        #endif
    }
    else
    {
        for (int i = world_rank - 1; i < N; i += world_size - 1)
        {        
            for (int j = 0; j < N; j++)
            {
                int aux[2], x, y;
                MPI_Recv(&aux, 2, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);                    
                x = aux[0];
                y = aux[1];
                C[i][j] = 0;
                for (int k = 0; k < N; k++)
                    C[i][j] += A[x][k] * B[y][k];
            }
        }

        for (int i = world_rank - 1; i < N;  i += world_size - 1)
            MPI_Send(C[i], N, MPI_INT, 0, 2, MPI_COMM_WORLD); 
    }

    // ____________________________________________

    MPI_Finalize();

    return 0;
}