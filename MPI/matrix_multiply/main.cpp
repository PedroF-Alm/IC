#include <stdio.h>
#include <time.h>
#include "./include/Mult.h"

#define N 1000
#define MAX 20
#define MIN 0

int main(int argc, char** argv) {
    
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

    srand(world_rank*time(NULL)); 

    if (world_rank == 0) {

        int **A = NULL, **B = NULL, **C = NULL;

        A = allocatemat(N);
        B = allocatemat(N);
        C = allocatemat(N);

        int min = MIN, max = MAX;

        if (argc > 2) {
            min = atoi(argv[1]);
            max = atoi(argv[2]);
        }
        
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = (int) ((double) rand() / RAND_MAX * (max - min + 1)) + min; 
                B[i][j] = (int) ((double) rand() / RAND_MAX * (max - min + 1)) + min;
            }
        }

        for (int i = 0; i < N; i++)
        {
            MPI_Send(A[i], N, MPI_INT, 1 + i % (world_size - 1), 1, MPI_COMM_WORLD);
            for (int j = 0; j < N; j++) {
                MPI_Send(B[j], N, MPI_INT, 1 + i % (world_size - 1), 2, MPI_COMM_WORLD);
            }
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                MPI_Recv(&C[i][j], 1, MPI_INT, 1 + i % (world_size - 1), 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("A[%d][%d] = %d\tB[%d][%d] = %d\tC[%d][%d] = %d\n", i, j, A[i][j], j, i, B[i][j], i, j, C[i][j]); 
            }
        }

        freemat(A, N);
        freemat(B, N);
        freemat(C, N);
        A = NULL;
        B = NULL;
        C = NULL;
    }
    else {

        int *A = NULL, *B = NULL, result;

        for (int i = world_rank; i < N; i += world_size - 1) {
            A = (int*) calloc(N, sizeof(int));
            MPI_Recv(A, N, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for (int j = 0; j < N; j++) {
                B = (int*) calloc(N, sizeof(int));
                MPI_Recv(B, N, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                result = multiply(A, B, N);
                MPI_Send(&result, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);

                if (B != NULL)
                    free(B), B = NULL;
            }

            if (A != NULL)
                free(A), A = NULL;
        }

    }

    // ____________________________________________

    MPI_Finalize();

    return 0;
}