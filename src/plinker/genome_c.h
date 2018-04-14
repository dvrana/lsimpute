
#include <stdio.h>

// struct fields are WIP, use interfacing functions below
typedef struct snp {
} snp_t;

typedef struct genome {
} genome_t;

// possibly may change
genome_t* g_fromfile(FILE* ped, FILE* map);

// number of individuals
int g_nsample(genome_t* g);

// number of SNPs
int g_nsnp(genome_t* g)

// Removes all SNPs from genome g not present in filt
void* g_filterby(genome_t* g, genome_t* filt);

// Filters out genomes not present among n given people in array ids
void* g_filterindiv(genome_t* g, int* ids, int n);

// Filters out all 
// Dylan's note to Cam, delete when done:
// Humans have 22 autosomes (present in everyone, two copies each) plus sex
// chromosomes and mitochondrial DNA.  That's complex.  For a first pass, remove
// everything not an autosome (chromosomes 1-22).
void* g_filterchrom(genome_t* g, int chromosome);

// What do these do?
snp_t *g_plookup(genome_t* g, int id);
snp_t *g_slookup(genome_t* g, int id);

// Recombinant distance (in centimorgans) between x and y
double s_pairdst(snp_t* x, snp_t* y);

