/* Interface for imputing samples from Li-Stephens probabilities
 */

#ifndef IMPUTE_H
#define IMPUTE_H

//#include <plinker/genome_c.h>

/* Returns the MLE for all ordered SNPs (according to ref) in sample.  It is
 * assumed that
 * a) SNPs in sample are a subset of ref
 * b) the nsample samples in ref are the same as the nsample columns in P
 * c) the nsnp SNPs in sample are the same as the nsnp rows in P
 */
snp_t* impute(float* P, genome_t* ref, genome_t* sample, int nsample);

#endif /* IMPUTE_H */
