/* Interface for Li-Stephens imputation.
 */

#include "../plinker/genome_c.h"
#include <stdlib.h>
#include <math.h>

#define EMISS(o1, o2, g) log(o1 == o2 ? (1 - g) : g)

// Takes the n log-scaled float values in A and returns their log-scaled sum
float logsum(float* A, int n) {
  return 0.0f; // TODO: this
}

// Adds two log-scaled probabilities
float logadd(float x, float y) {
  return 0.0f; // TODO: this
}

// Calculates ln(1.0f - x), where x is a log-scaled probability
float logsub1(float x) {
  return 0.0f; // TODO: this
}


/* Forward algorithm
 * Note that the memory it returns (same form as ls) is heap-allocated and so
 * must be freed.  fw[i][j] is the probability that we are in the jth state
 * given SNPs [0,i].
 * Note also that the probabilities it returns are ln-scaled.
 */
float* forward(genome_t sample, std::string id, genome_t ref, float g, float theta) {
  int n_ref = g_nsample(ref);
  int n_snp = g_nsnp(ref);

  // Initialize memory
  float* fw = (float*)malloc(sizeof(float) * n_snp * n_ref);
  
  // Get snp arrays
  snp_t* s = g_plookup(sample,id);
  snp_t** S = new snp_t*[n_ref];
  int i = 0;
  for (auto id : *ref) {
    S[i] = g_plookup(ref,id.first);
    i++;
  }

  // Initialize the first row
  float c = log(1.0f / ((float)n_ref)); // probability of jumping to given ref
  for (int i = 0; i < n_ref; i++) fw[i] = c + EMISS(s[0], S[i][0], g);

  // For each iteration
  for (int i = 1; i < n_snp; i++) {
    // Precompute jump probability
    float J = logsum(&(fw[(i-1)* n_snp]),n_ref);
    J = J + log((1.0f - exp(-1 * theta * g_rec_dist(ref, i-1))));
    float nJ = logsub1(J);
    
    // Calculate values
    for (int j = 0; j < n_ref; j++) {
      float alpha = logadd((fw[((i-1)*n_snp)+j] + nJ),(J + c));
      fw[i * n_snp + j] = alpha + EMISS(s[i],S[j][i],g);
    }
  }
  return fw;
}

/* Backward algorithm
 * Note that the memory it returns (same form as ls) is heap-allocated and so
 * must be freed. bw[i][j] is the probability that we observe SNPs [i,n]
 * given that we are in state j at time i.
 * Note also that the probabilities it returns are ln-scaled.
 */
float* backward(genome_t sample, std::string id, genome_t ref, float g, float theta) {
  int n_ref = g_nsample(ref);
  int n_snp = g_nsnp(ref);

  // Initialize memory
  float* bw = (float*)malloc(sizeof(float) * n_snp * n_ref);
  
  // Get snp arrays
  snp_t* s = g_plookup(sample,id);
  snp_t** S = new snp_t*[n_ref];
  int i = 0;
  for (auto id : *ref) {
    S[i] = g_plookup(ref,id.first);
    i++;
  }

  // Initialize the last row
  float c = log(1.0f / ((float)n_ref)); // probability of jumping to given ref
  for (int i = (n_snp - 1) * n_ref; i < (n_snp * n_ref); i++)
    bw[i] = EMISS(s[n_snp - 1],S[i][n_snp - 1],g);

  // For each iteration
  for (int i = n_snp - 2; i >= 0; i--) {
    // Precompute reverse jump probability
    float J = logsum(&(bw[(i+1)* n_snp]),n_ref);
    J = J + log((1.0f - exp(-1 * theta * g_rec_dist(ref, i))));
    float nJ = logsub1(J);

    // Calculate values
    for (int j = 0; j < n_ref; j++) {
      bw[i * n_ref + j] =  
        logadd(J + c, nJ + bw[(i+1) * n_ref + j]) + EMISS(s[i],S[j][i],g);
    }
  }
  return bw;
}

/* Returns smoothed Li-Stephens probabilities as a two-dimensional, 
 * heap-allocated array A[s][n], where s is the number of SNPs and n the number
 * of reference genomes, and A[i][j] is the natural log of the probability that
 * the ancestor is reference genome j at SNP i.
 *
 * Assumes that SNPs present in reference and sample are the same.
 *
 * Remember to call free on the result.
 * 
 * Probabilities are for genome id from sample, using ref as a reference panel.
 * g  - Garble rate- probability that a test haplotype doesn't line up with
 *   reference from which it comes.
 * theta - Recombination rate constant, s.t. jump probability is 
 *   1-e^{-theta d}, where d is distance in centimorgans
 */
float* ls(genome_t sample, std::string id, genome_t ref, float g, float theta) {
  // Forward pass
  float* fw = forward(sample, id, ref, g, theta);
  
  // Backward pass
  float* bw = backward(sample, id, ref, g, theta);

  // Smoothing pass
  // TODO: this

  return NULL; // TODO: this
}