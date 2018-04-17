
#include <stdio.h>
#include <unordered_map>
#include <memory>

// struct fields are WIP, use interfacing functions below
typedef struct snp {
} snp_t;

struct genome {
    int nsnp;
    int nsample;
    // TODO: use more than individual id
    std::unordered_map<int, snp_t*> samples;
};

// NOTE TO DYLAN (delete):
//   I'm choosing to use shared_ptrs to avoid the hassle of having to manually
//   manage the memory. Let me know if this causes issues -- it shouldn't,
//   beyond possibly needing to change function headers.
typedef std::shared_ptr<struct genome> genome_t;

genome_t g_empty();
genome_t g_fromfile(FILE* ped, FILE* map);

// number of individuals
int g_nsample(genome_t g);

// number of SNPs
int g_nsnp(genome_t g);

// Removes all SNPs from genome g not present in filt
void g_filterby(genome_t g, genome_t filt);

// Filters out genomes not present among n given people in array ids
void g_filterindiv(genome_t g, int* ids, int n);

// Humans have 22 autosomes (present in everyone, two copies each) plus sex
// chromosomes and mitochondrial DNA. That's complex. For a first pass, remove
// everything not an autosome (chromosomes 1-22).
void g_filterchrom(genome_t g, int chromosome);

// Lookup SNP list by person
snp_t* g_plookup(genome_t g, int pid);

// Lookup SNP by person and SNP identifier
snp_t g_ilookup(genome_t g, int pid, int sid);

// Genetic distance (in centimorgans) of s
double s_dst(snp_t s);

