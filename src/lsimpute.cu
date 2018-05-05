
#include <cuda.h>
#include "lsimpute.h"

#define EMISS(o1, o2, g) log(o1 == o2 ? (1 - g) : g)

#define BLOCKMAX 1024 // TODO: set to 1024 if compute capability >= 2.0

/* Calculates log(exp(x) + exp(y))
 */
__device__ float d_logadd(float x, float y) {
  if (y > x) { return d_logadd(y,x); }
  return x + logf(1+expf(y-x));
}

/* Calculates log(1.0f - exp(x))
 */
__device__ float d_logsub1(float x) {
  return logf(1.0f - expf(x));
}

__device__ bool isPow2(int n) {
  return (n != 0) && (n & (n-1)) == 0;
}

/* Sums the n floating point values in array A (in no particular order)
 * that are in natural log space
 * Basically, find log(sum(p1 ... pn)) given log(p1) .. log(pn)
 */
__device__ float row_logsum(float* A, int n, float* scratch) {
  int tid = threadIdx.x;
  int i = blockIdx.x*blockDim.x + threadIdx.x;

  int elts = n / BLOCKMAX + (tid < n % BLOCKMAX);

  // With sample sizes being single-digit multiples of BLOCKMAX (at least at
  // the moment), we can get away with doing this linearly. As sample sizes
  // get larger, we'll need to start looking into recursive kernel invocation.
  scratch[tid] = A[i];
  for (int j = 1 ; j < elts ; j += 1) {
    scratch[tid] = d_logadd(scratch[tid], A[i+j]);
  }
  __syncthreads();

  for (int s = 1 ; s < blockDim.x ; s <<= 1) {
    int index = 2*s*tid;
    if (index < blockDim.x) {
      scratch[tid] = d_logadd(scratch[index], scratch[index+s]);
    }
  }
  __syncthreads();

  return scratch[0];
}

/* Normalizes the n floating point values in the array starting at A
 * (sets A[i] = log(exp(A[i]) / exp(reduce_logsum(A,n))))
 */
__device__ void d_logrownorm(float* A, int n, float* scratch) {
  float x = row_logsum(A, n, scratch);
  for (int i = 0; i < n; i++) A[i] -= x;
  return;
}

__device__ void fwKernel(uint8_t* refs, uint8_t* sample, float* dists,
    float* fw, float g, float theta, int nsnp, int nsample, float* scratch) {
  // Initialize first row
  float c = 1.0f / ((float)nsample);
  for (int i = threadIdx.x; i < nsample; i += blockDim.x) fw[i] = c;
  __syncthreads();

  // For each SNP (going forward)
  for (int k = 1; k < nsnp; k++) {
    int K = k * nsample;
    // Precompute jump probability
    int J = row_logsum(&(fw[K]), nsample, scratch);
    J = J + d_logsub1(-1.0f * theta * dists[k]);
    float nJ = d_logsub1(J);

    // Calculate values
    for (int i = threadIdx.x; i < nsample; i += blockDim.x) {
      float alpha = d_logadd(fw[K - nsample + i] + nJ, J + c);
      fw[K + i] = alpha + EMISS(sample[k], refs[K+i],g);
    }
  }
  return;
}

__device__ void bwKernel(uint8_t* refs, uint8_t* sample, float* dists,
    float* bw, float g, float theta, int nsnp, int nsample, float* scratch) {
  // Initialize last row
  float c = 1.0f / ((float)nsample);
  for (int i = threadIdx.x; i < nsample; i += blockDim.x) {
    bw[(nsample * (nsnp - 1)) + i] =
      c + EMISS(refs[nsample * (nsnp - 1) + i], sample[nsnp - 1], g);
  }

  // For each SNP (going backward)
  for (int k = nsnp - 2; k >= 0; k--) {
    int K = k * nsample;
    // Precompute jump probability
    float J = row_logsum(&(bw[K+nsample]),nsample, scratch);
    J = J + d_logsub1(-1.0f * theta * dists[k]);
    float nJ = d_logsub1(J);

    // Calculate values
    for (int i = threadIdx.x; i < nsample; i += blockDim.x) {
      float alpha = d_logadd(J + c, nJ + bw[K + i + nsample]);
      bw[K + i] = alpha + EMISS(sample[k], refs[K+i], g);
    }
  }
  return;
}

/* Returns its smoothed answer in fw
 */
__device__ void smoothKernel(
    float* fw, float* bw, int nsnp, int nsample, float* scratch
) {
  for (int i = 0; i < nsnp; i++) {
    int I = i * nsample;
    for (int j = 0; j < nsample; j++) {
      fw[I + j] = fw[I + j] + bw[I + j];
    }
    d_logrownorm(&(fw[I]), nsample, scratch);
  }
  return;
}

__global__ void computeKernel(uint8_t* refs, uint8_t* sample, float* dists,
    float* fw, float* bw, float g, float theta, int nsnp, int nsample) {
  extern __shared__ float scratch[];
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

  // Run the kernel
  computeKernel<<<1, BLOCKMAX, nsample*sizeof(float)>>>
    (d_refs, d_sample, d_dists, d_fw,
      d_bw, g, theta, nsnp, nsample);
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

