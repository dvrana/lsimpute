/* Interface for Li-Stephens imputation.
 */

#ifndef LS_H
#define LS_H

//#include <plinker/genome_c.h>

/* Returns smoothed Li-Stephens probabilities as a two-dimensional, 
 * heap-allocated array A[s][n], where s is the number of SNPs and n the number
 * of reference genomes, and A[i][j] is the natural log of the probability that
 * the ancestor is reference genome j at SNP i.
 * 
 * Probabilities are for genome id from sample, using ref as a reference panel.
 * g  - Garble rate- probability that a test haplotype doesn't line up with
 *   reference from which it comes.
 * theta - Recombination rate constant, s.t. jump probability is 
 *   1-e^{-theta d}, where d is distance in centimorgans
 */
float* ls(genome_t sample, std::string id, genome_t ref, float g, float theta); 

#endif /* LS_H */
