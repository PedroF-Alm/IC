#define N 10000000
#include <stdio.h>

__global__ void vector_add(float *out, float *a, float *b, int n) {
    int index = threadIdx.x;
    int stride = blockDim.x;
    for(int i = index; i < n; i += stride){
        out[i] = a[i] + b[i];
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
    vector_add<<<1, 256>>>(d_out, d_a, d_b, N);

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