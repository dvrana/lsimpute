
// HEADER FILE FOR GPU CODE

#ifndef LSIMPUTE_CU_H
#define LSIMPUTE_CU_H

#include <cstdint>
#include "plinker/genome_c.h"

class lsimputer {
    // XXX: coding style: it would be more idiomatic to make these private,
    // and to follow proper OO design patterns. However, I prefer to use
    // these as glorified structs, so for now we'll make everything public.
public:
    int nsnp;
    int nsample;
    // SNPs are stored in SNP-major order. The ordering of samples is based on
    // the internal ordering of the input genome.
    // For now, SNPs are represented using full 8-bit valued bitvectors. If
    // memory is a concern, we can fix this.
    uint8_t* ref;
    double* dists;

    // Constants for model
    float g;
    float theta;

    lsimputer(genome_t g);

    ~lsimputer();

    // This function should probably be defined in the CUDA files, by defining
    // the function lsimpute::compute in lsimpute.cu (or equivalents)
    // Returns a malloc'ed float array
    float* compute(uint8_t* snps);
};

#endif

