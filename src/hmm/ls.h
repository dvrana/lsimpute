/* Interface for Li-Stephens imputation.
 */

#include <plinker/genome_c.h>

/* Returns smoothed Li-Stephens probabilities as a two-dimensional, 
 * heap-allocated array A[s][n], where s is the number of SNPs and n the number
 * of reference genomes, and A[i][j] is the natural log of the probability that
 * the ancestor is reference genome j at SNP i.
 * 
 * Probabilities are for genome id from sample, using ref as a reference panel.
 * g  - TODO: WHAT IS
 * ck - TODO: WHAT IS
 * theta - Recombination rate constant
 */
float* ls(genome_t* sample, int id, genome_t* ref, float g, float ck, 
    float theta);
