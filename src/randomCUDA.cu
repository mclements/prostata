#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <stdio.h>

// 1. Initialization Kernel
// Sets up a distinct, independent stream for each thread using MRG32k3a
__global__ void setup_mrg_kernel(curandStateMRG32k3a *state, unsigned long long seed) {
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    
    // Each thread gets the same seed, but a different sequence/stream number (id)
    // parameter 3 (offset) is set to 0
    curand_init(seed, id, 0, &state[id]);
}

// 2. Generation Kernel
__global__ void generate_mrg_kernel(curandStateMRG32k3a *state, double *output, int values_per_thread) {
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    
    // Cache the state in local registers for maximum performance
    curandStateMRG32k3a localState = state[id];
    
    int base_index = id * values_per_thread;
    
    for (int i = 0; i < values_per_thread; i++) {
        // cuRAND outputs (0.0, 1.0] for double precision
        double raw_rand = curand_uniform_double(&localState);
        
        // Transform to [0.0, 1.0) to match dSFMT's close_open property
        output[base_index + i] = 1.0 - raw_rand;
    }
    
    // Save the updated state back to global memory for future kernel calls
    state[id] = localState;
}

void test_mrg32k3a_random_generation() {
    // Thread configuration
    int threadsPerBlock = 256;
    int blocksPerGrid = 100;
    int totalThreads = threadsPerBlock * blocksPerGrid;
    int valuesPerThread = 10;
    int totalValues = totalThreads * valuesPerThread;

    // Allocate Device memory for cuRAND states
    curandStateMRG32k3a *d_states;
    cudaMalloc((void**)&d_states, totalThreads * sizeof(curandStateMRG32k3a));

    // Allocate Device memory for output data
    double *d_output;
    cudaMalloc((void**)&d_output, totalValues * sizeof(double));

    // 1. Initialize the MRG32k3a states on the GPU
    unsigned long long seed = 1234ULL;
    setup_mrg_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_states, seed);
    cudaDeviceSynchronize();

    // 2. Generate the [0.0, 1.0) random doubles
    generate_mrg_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_states, d_output, valuesPerThread);
    cudaDeviceSynchronize();

    // Allocate Host memory to verify results
    double *h_output = (double*)malloc(totalValues * sizeof(double));
    cudaMemcpy(h_output, d_output, totalValues * sizeof(double), cudaMemcpyDeviceToHost);

    // Print a few sample values to verify bounds
    printf("Sample outputs (should be inside [0.0, 1.0)):\n");
    for (int i = 0; i < 5; i++) {
        printf("Value[%d]: %f\n", i, h_output[i]);
    }

    // Cleanup
    cudaFree(d_states);
    cudaFree(d_output);
    free(h_output);

}
