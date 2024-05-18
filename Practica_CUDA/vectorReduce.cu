/*
 * This example shows how to compute the reduction of the elements of a vector.
 *
 * Also, it shows how to measure the performance of block of threads of a kernel
 * accurately. Blocks are executed in parallel and out of order. Since there's no
 * synchronization mechanism between blocks, we measure the clock once for each block.
 *
 */

// Includes C
#include <stdio.h>
#include <assert.h>

// CUDA Runtime
#include <cuda_runtime.h>

// Includes Helper Functions
#include <helper_cuda.h>
#include <helper_functions.h>

/*
 * vectorReduce
 *
 * This kernel computes a standard parallel reduction and evaluates the
 * time it takes to do that for each block. The timing results are stored in device memory.
 * 
 */
__global__ void time_and_reduce(float *vector_d, float *reduce_d, clock_t *times_d, int n)
{
	extern __shared__ float sdata[];
	
	// local thread ID (in block)
	int tidb = threadIdx.x;
	
    // global thread (ID in grid)
	int tidg = blockIdx.x * blockDim.x + tidb;
	
	// record the initial time for each block
	if (tidb == 0) {
		times_d[blockIdx.x] = clock();
	}
	
	// move data from global to shared memory
	sdata[tidb] = (tidg < n) ? vector_d[tidg] : 0;
	__syncthreads();
	
	// perform reduction in shared memory
	for (int s = blockDim.x/2; s > 0; s >>= 1) {
		if (tidb < s) {
			sdata[tidb] += sdata[tidb + s];
		}
		__syncthreads();
	}
	
	// write result for this block to global memory
	if (tidb == 0) {
		reduce_d[blockIdx.x] = sdata[0];
		times_d[blockIdx.x] = clock() - times_d[blockIdx.x];
	}
}

/*
 * Host main routine
 *
 */
int main(int argc, char **argv)
{
	// default parameter values
	int n = 1024, bsx = 32;
	
	// process command line arguments
	if (checkCmdLineFlag(argc, (const char **)argv, "n")) {
		n = getCmdLineArgumentInt(argc, (const char **)argv, (const char *) "n");
	}
	if (checkCmdLineFlag(argc, (const char **)argv, "bsx")) {
		bsx = getCmdLineArgumentInt(argc, (const char **) argv, (const char *) "bsx");
	}
	size_t nBytes = n * sizeof(float);
	
	clock_t *clocks_h = NULL;
	clock_t *clocks_d = NULL;
	
	float elapsed_time = .0;
	float *vector_h, *reduce_h;	// host data
    float *vector_d, *reduce_d;	// device data
	
	// set the GPU to use
	int dev = 0;
	cudaSetDevice(dev);
	
	// total number of thread blocks
	int nblocks = (n + bsx - 1) / bsx;

	// set kernel launch configuration
    dim3 grid(nblocks);
    dim3 block(bsx);
	
    // allocate host memory
    vector_h = (float *) malloc(nBytes);
    clocks_h = (clock_t *) malloc(nblocks * sizeof(clock_t));
    reduce_h = (float *) malloc(nblocks * sizeof(float));
	
	float acum = .0;
	// initialize host memory
    for(int i = 0; i < n; i++) {
        vector_h[i] = (float) 1;
		acum += 1.0;
	}
	
    // allocate device memory
    cudaMalloc((void **) &vector_d, nBytes);
    cudaMalloc((void **) &reduce_d, nblocks * sizeof(float));
	cudaMalloc((void **) &clocks_d, nblocks * sizeof(clock_t));
	
	// create cuda events
    cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	// insert stream 0 in start event
	cudaEventRecord(start, 0);
	
    // copy data from host memory to device memory
    cudaMemcpy(vector_d, vector_h, nBytes, cudaMemcpyHostToDevice);
    
    // execute the kernel 
    printf("---> Running configuration: grid of %d blocks of %d threads (TOTAL: %d threads)\n", nblocks, bsx, nblocks * bsx);
    time_and_reduce<<<grid, block, bsx * sizeof(float)>>>(vector_d, reduce_d, clocks_d, n);

    // copy data from device memory to host memory
	cudaMemcpy(reduce_h, reduce_d, nblocks * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(clocks_h, clocks_d, nblocks * sizeof(clock_t), cudaMemcpyDeviceToHost);
	
	// insert stream 0 in stop event
	cudaEventRecord(stop, 0);

    // using events to calculate the execution time        
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed_time, start, stop);
    printf("---> Time spent executing by the GPU: %.2f\n", elapsed_time);
	
	long double avgElapsedClocks = 0;
    for (int i = 0; i < nblocks; i++) {
		avgElapsedClocks += (long double) clocks_h[i];
    }
    avgElapsedClocks = avgElapsedClocks / nblocks;
    printf("Average Clocks/Block = %Lf\n", avgElapsedClocks);

	// check the output for correctness
	float result = 0.0;
	for(int i = 0; i < nblocks; i++) { 
		result += reduce_h[i];
	}
	assert(result == (float) acum);
	
	// destroy events
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

    // free host memory
    free(vector_h);
	free(reduce_h);
	free(clocks_h);
	
	// free device memory
    cudaFree(vector_d);
    cudaFree(reduce_d);
	cudaFree(clocks_d);
	
    printf("\nTest PASSED\n");
	exit(EXIT_SUCCESS);
}
