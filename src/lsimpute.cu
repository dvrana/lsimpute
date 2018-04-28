
#include "lsimpute.h"

#define BLOCKMAX 512 // TODO: set to 1024 if compute capability >= 2.0

__global__ void computeKernel(uint8_t* refs, uint8_t* sample, float* dists,
    float* fw, float* bw, float g, float theta) {
  return; // TODO: this
}

float* lsimputer::compute(uint8_t* snps) {
  // Allocate space for refs, sample, distances, and return values
  uint8_t* d_refs;
  uint8_t* d_sample;
  float* d_dists;
  float* d_fw;
  float* d_bw;
  cudaMalloc((void **)&d_refs, sizeof(uint8_t) * nsnp * nsample);
  cudaMalloc((void **)&d_sample, sizeof(uint8_t) * nsnp);
  cudaMalloc((void **)&d_dists, sizeof(float) * nsnp);
  cudaMalloc((void **)&d_fw, sizeof(float) * nsnp * nsample);
  cudaMalloc((void **)&d_bw, sizeof(float) * nsnp * nsample);

  // Transfer over data
  cudaMemcpy(d_refs, ref, sizeof(uint8_t) * nsnp * nsample,
      cudaMemcpyHostToDevice);
  cudaMemcpy(d_sample, snps, sizeof(uint8_t) * nsnp,
      cudaMemcpyHostToDevice);
  cudaMemcpy(d_dists, dists, sizeof(float) * nsnp,
      cudaMemcpyHostToDevice);

  // Run the kernel
  computeKernel<<<1, BLOCKMAX>>>(d_refs, d_sample, d_dists, d_fw, d_bw, g, theta);
  cudaDeviceSynchronize();

  // Transfer data off the device
  float* res = (float*)malloc(sizeof(float) * nsnp * nsample);
  cudaMemcpy(res, d_fw, sizeof(uint8_t) * nsnp * nsample,
      cudaMemcpyDeviceToHost);

  // Free device memory
  cudaFree(d_refs);
  cudaFree(d_sample);
  cudaFree(d_dists);
  cudaFree(d_fw);
  cudaFree(d_bw);

  return res;
}

