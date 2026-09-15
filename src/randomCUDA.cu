#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <stdio.h>
#include <stdexcept>
#include <sstream>

/***********************************************************************************************************************************
 *
 * This macro checks return value of the CUDA runtime call and
 * throws an exception if the call failed.
 *
 * @param value cuda function call that returns cudaError_t value
 */
#define gpuErrchk(ans) \
    { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char* file, int line) {
    if (code != cudaSuccess) {
        std::ostringstream errorStream;
        errorStream << "Error " << cudaGetErrorString(code) << " at line " << line << " in file " << file;
        throw std::runtime_error(errorStream.str());
    }
}

// 1. Initialization Kernel
// Sets up a distinct, independent stream for each thread using MRG32k3a
__global__ void initRNG(curandStateMRG32k3a *state, unsigned long long seed, int N) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid >= N) return;
    
    // Each thread gets the same seed, but a different sequence/stream number (id)
    // parameter 3 (offset) is set to 0
    curand_init(seed, tid, 0, &state[tid]);
}

// 2. Generation Kernel
__global__ void generate_mrg_kernel(curandStateMRG32k3a *state, double *output, int values_per_thread, int N) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid * values_per_thread + values_per_thread - 1 >= N) return;
    
    // Cache the state in local registers for maximum performance
    curandStateMRG32k3a localState = state[tid];
    
    for (int i = 0; i < values_per_thread; i++) {
        // cuRAND outputs (0.0, 1.0] for double precision
        double raw_rand = curand_uniform_double(&localState);
        
        // Transform to [0.0, 1.0) to match dSFMT's close_open property
        output[tid * values_per_thread + i] = 1.0 - raw_rand;
    }
    
    // Save the updated state back to global memory for future kernel calls
    state[tid] = localState;
}

cudaDeviceProp initDevice(int device_id) {
    cudaDeviceProp prop;
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0) {
        fprintf(stderr, "No CUDA devices found.\n");
        exit(EXIT_FAILURE);
    }
    if (cudaSetDevice(device_id) != cudaSuccess) {
        cudaDeviceReset();
        fprintf(stderr, "Failed to set CUDA device %d\n", device_id);
        exit(EXIT_FAILURE);
    }
    if (cudaGetDeviceProperties(&prop, device_id) != cudaSuccess) {
        cudaDeviceReset();
        fprintf(stderr, "Failed to get properties for CUDA device %d\n", device_id);
        exit(EXIT_FAILURE);
    }
    return prop;
}

void test_mrg32k3a_random_generation(int N, double *h_output) {
    
    cudaDeviceProp prop = initDevice(0);

    dim3 block, grid;
    block.x = 64;
    grid.x  = (N + block.x - 1) / block.x;
    
    // Aim to launch around ten or more times as many blocks as there
    // are multiprocessors on the target device.
    unsigned int blocksPerSM = 10;
    unsigned int numSMs      = prop.multiProcessorCount;
    
    while (grid.x > 2 * blocksPerSM * numSMs) {
        grid.x >>= 1;
    }
    int totalThreads = grid.x * block.x;
    int valuesPerThread = (N + totalThreads - 1) / totalThreads;

    // Allocate Device memory for cuRAND states
    curandStateMRG32k3a* d_rngStates;
    gpuErrchk(cudaMalloc((void**)&d_rngStates, N * sizeof(curandStateMRG32k3a)));
    unsigned long long seed = 1234ULL;
    initRNG<<<grid, block>>>(d_rngStates, seed, N);
    gpuErrchk(cudaDeviceSynchronize());

    // 2. Generate the [0.0, 1.0) random doubles
    double* d_output;
    gpuErrchk(cudaMalloc((void**)&d_output, N * sizeof(double)));
    generate_mrg_kernel<<<grid, block>>>(d_rngStates, d_output, valuesPerThread, N);
    gpuErrchk(cudaDeviceSynchronize());

    // Allocate Host memory to verify results
    gpuErrchk(cudaMemcpy(h_output, d_output, N * sizeof(double), cudaMemcpyDeviceToHost));

    cudaFree(d_rngStates);
    cudaFree(d_output);

}
