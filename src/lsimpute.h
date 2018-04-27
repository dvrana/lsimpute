
// HEADER FILE FOR GPU CODE

#ifndef LSIMPUTE_CU_H
#define LSIMPUTE_CU_H

#include <cstdint>

struct lsimputer {
    int nsnp;
    // For now, SNPs are represented using full 8-bit valued bitvectors. If
    // memory is a concern, we can fix this.
    uint8_t* snps;
    double* dists;

    // Constants for model
    float g;
    float theta;

    lsimputer(int nsnp_) {
        nsnp = nsnp_;
        snps = new uint8_t[nsnp];
        dists = new double[nsnp-1];
    }

    ~lsimputer() {
        delete snps;
        delete dists;
    }

    // This function should probably be defined in the CUDA files, by defining
    // the function lsimpute::compute in lsimpute.cu (or equivalents)
    // Returns a malloc'ed float array
    float* compute();
};

#endif

