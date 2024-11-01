#define N 10000000
#include <stdio.h>

__global__ void vector_add(float *out, float *a, float *b, int n) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Handling arbitrary vector size
    if (tid < n){
        out[tid] = a[tid] + b[tid];
    }
}

int main(){
    float *h_a, *h_b, *h_out, *d_a, *d_b, *d_out; 

    // Allocate memory (host)
    h_a   = (float*)malloc(sizeof(float) * N);
    h_b   = (float*)malloc(sizeof(float) * N);
    h_out = (float*)malloc(sizeof(float) * N);

    // Allocate memory (device)
    cudaMalloc((void**)&d_a,   sizeof(float) * N);
    cudaMalloc((void**)&d_b,   sizeof(float) * N);
    cudaMalloc((void**)&d_out, sizeof(float) * N);

    // Initialize arrays
    for(int i = 0; i < N; i++){
        h_a[i] = 1.0f; h_b[i] = 2.0f;
    }

    // Transfer data from host to device memory
    cudaMemcpy(d_a,   h_a,   sizeof(float) * N, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b,   h_b,   sizeof(float) * N, cudaMemcpyHostToDevice);
    cudaMemcpy(d_out, h_out, sizeof(float) * N, cudaMemcpyHostToDevice);

    // Main function
    int block_size = 256;
    int grid_size = ((N + block_size) / block_size);
    vector_add<<<grid_size,block_size>>>(d_out, d_a, d_b, N);

    // Transfer data from device to host memory
    cudaMemcpy(h_out, d_out, sizeof(float) * N, cudaMemcpyDeviceToHost);

    printf("%f\n", h_out[0]);

    // Cleanup after kernel execution
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_out);
    free(h_a);
    free(h_b);
    free(h_out);

    // out array must be transfered back and forth
}