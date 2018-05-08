
#include <cuda.h>
#include <math_constants.h>
#include "lsimpute.h"

#define EMISS(o1, o2, g) (o1 == o2 ? (1 - g) : g)

#define BLOCKMAX 512 // TODO: set to 1024 if compute capability >= 2.0
#define WARP_SIZE 32

__device__ bool isPow2(int n) {
  return (n != 0) && (n & (n-1)) == 0;
}

__device__ __host__
int npow2(int n) {
  n -= 1;
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
  return n+1;
}

/* Sums the n floating point values in array A (in no particular order)
 * that are in natural log space
 * Basically, find log(sum(p1 ... pn)) given log(p1) .. log(pn)
 */
__device__ float row_sum(float* A, int n, float* scratch) {
  n = npow2(n);
  // Possible performance speedup: half our threads are idle on the first loop!
  int tid = threadIdx.x;
  int nthread = blockDim.x;
  int fakeelts = n / nthread;

  int i = blockIdx.x*blockDim.x + (fakeelts * tid) + min(n % nthread, tid);
  int elts = fakeelts + (tid < n % nthread);

  // With sample sizes being single-digit multiples of BLOCKMAX (at least at
  // the moment), we can get away with doing this linearly. As sample sizes
  // get larger, we'll need to start looking into recursive kernel invocation.
  scratch[tid] = A[i];
  for (int j = 1 ; j < elts && i+j < n ; j += 1) {
    scratch[tid] = scratch[tid] + A[i+j];
  }
  __syncthreads();

  int big = n > WARP_SIZE;

  for (int s = nthread/2 ; s > WARP_SIZE*big ; s >>= 1) {
    if (tid < s && tid+s < n) {
      scratch[tid] = scratch[tid] + scratch[tid+s];
    }
    __syncthreads();
  }

  if (tid < WARP_SIZE && big) {
    #pragma unroll
    for (int j = WARP_SIZE ; j > 0 ; j >>= 1) {
      if (tid+j < n)
        scratch[tid] = scratch[tid] + scratch[tid+j];
    }
  }

  __syncthreads();
  return scratch[0];
}

/* Normalizes the n floating point values in the array starting at A
 * (sets A[i] = log(exp(A[i]) / exp(reduce_logsum(A,n))))
 */
__device__ void d_rownorm(float* A, int n, float* scratch) {
  float x = row_sum(A, n, scratch);
  for (int i = 0; i < n; i++) A[i] /= x;
  return;
}

__device__ void fwKernel(uint8_t* refs, uint8_t* sample, float* dists,
    float* fw, float g, float theta, int nsnp, int nsample, float* scratch) {
  // Initialize first row
  float c = 1.0f / ((float)nsample);
  for (int i = threadIdx.x; i < nsample; i += blockDim.x)
    fw[i] = EMISS(sample[0], refs[i],g);
  __syncthreads();

  // For each SNP (going forward)
  for (int k = 1; k < nsnp; k++) {
    int K = k * nsample;
    // Precompute jump probability
    float x = row_sum(&(fw[K-nsample]), nsample, scratch);
    float nJ = expf(-1.0f * theta * dists[k-1]);
    float J = 1.0 - nJ;

    // Calculate values
    for (int i = threadIdx.x; i < nsample; i += blockDim.x) {
      fw[K - nsample + i] /= x;
      float alpha = fw[K - nsample + i] * nJ + J * c;
      fw[K + i] = alpha * EMISS(sample[k], refs[K+i],g);
    }
    __syncthreads();
  }
  return;
}

__device__ void bwKernel(uint8_t* refs, uint8_t* sample, float* dists,
    float* bw, float g, float theta, int nsnp, int nsample, float* scratch) {
  // Initialize last row
  float c = 1.0f / ((float)nsample);
  for (int i = threadIdx.x; i < nsample; i += blockDim.x) {
    bw[(nsample * (nsnp - 1)) + i] =
      EMISS(refs[nsample * (nsnp - 1) + i], sample[nsnp - 1], g);
  }

  // For each SNP (going backward)
  for (int k = nsnp - 2; k >= 0; k--) {
    int K = k * nsample;
    // Precompute jump probability
    float x = row_sum(&(bw[K+nsample]), nsample, scratch);
    float nJ = expf(-1.0f * theta * dists[k]);
    float J = 1.0f - nJ;

    // Calculate values
    for (int i = threadIdx.x; i < nsample; i += blockDim.x) {
      bw[K + nsample + i] /= x;
      float alpha = (J * c) + (nJ * bw[K + i + nsample]);
      bw[K + i] = alpha * EMISS(sample[k], refs[K+i], g);
    }
    __syncthreads();
  }
  return;
}

/* Returns its smoothed answer in fw
 */
__device__ void smoothKernel(
    float* fw, float* bw, int nsnp, int nsample, float* scratch) {
  for (int i = 0; i < nsnp-1; i++) {
    int I = i * nsample;
    for (int j = threadIdx.x; j < nsample; j += blockDim.x) {
      fw[I + j] = fw[I + j] * bw[I + nsample + j];
    }
    __syncthreads();
    d_rownorm(fw + I, nsample, scratch);
  }
  d_rownorm(fw + (nsample * (nsnp-1)), nsample, scratch);
  return;
}

__global__ void computeKernel(uint8_t* refs, uint8_t* sample, float* dists,
    float* fw, float* bw, float g, float theta, int nsnp, int nsample,
    int nscratch) {
  extern __shared__ float scratch[];

  __syncthreads();

  // Forward step
  fwKernel(refs, sample, dists, fw, g, theta, nsnp, nsample, scratch);

  // Backward step
  bwKernel(refs, sample, dists, bw, g, theta, nsnp, nsample, scratch);

  // Smoothing step
  smoothKernel(fw, bw, nsnp, nsample, scratch);
  return;
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

  int nthread = npow2(min(BLOCKMAX, max(nsample, 32)));
  int nscratch = max(nthread, npow2(nsample));

  // Run the kernel
  computeKernel<<<1, nthread, nscratch*sizeof(float)>>>
    (d_refs, d_sample, d_dists, d_fw,
      d_bw, g, theta, nsnp, nsample, nscratch);
  cudaDeviceSynchronize();

  // Transfer data off the device
  float* res = (float*)malloc(sizeof(float) * nsnp * nsample);
  cudaMemcpy(res, d_fw, sizeof(float) * nsnp * nsample,
      cudaMemcpyDeviceToHost);

  // Free device memory
  cudaFree(d_refs);
  cudaFree(d_sample);
  cudaFree(d_dists);
  cudaFree(d_fw);
  cudaFree(d_bw);

  return res;
}

