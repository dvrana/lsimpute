
#include <stdio.h>

// struct fields are WIP, use interfacing functions below
typedef struct snp {
} snp_t;

typedef struct genome {
} genome_t;

// possibly may change
genome_t *g_fromfile(FILE *ped, FILE *map);

void *g_filterby(genome_t *g, genome_t *filt);

snp_t *g_plookup(genome_t *g, int id);
snp_t *g_slookup(genome_t *g, int id);

// possibly might want to use a different type
size_t s_absdst(snp_t *);
size_t s_pairdst(snp_t *, snp_t *);

