
#ifndef GENOME_C_H
#define GENOME_C_H

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <bitset>
#include <iterator>

enum allele { A, C, G, T };

typedef allele snp_t;

struct snpmeta {
    int ind; // index in mapfile ordering
    std::string id;
    int chnum;
    double gdist;
    int pos;

    bool operator < (const snpmeta& s) const {
        return pos < s.pos;
    }
};

struct snpmap {
    int nsnp;
    // maps index in mapfile to index in bp ordering
    std::shared_ptr<int> ids;
    // maps snp id string to index in bp ordering
    std::shared_ptr<std::string> sids;
    // ordered in bp order
    std::shared_ptr<std::vector<struct snpmeta>> data;

    int* id_arr() {
        return ids.get();
    }
};

struct genome {
    int nsample;
    struct snpmap map;
    // maps familyid_indid to sample
    std::map<std::string, std::shared_ptr<snp_t>> samples;

    typedef std::map<std::string, std::shared_ptr<snp_t>>::iterator iter;
    iter begin() { return samples.begin(); }
    iter end() { return samples.end(); }
};

typedef std::shared_ptr<struct genome> genome_t;

struct genomeErr : public std::exception {
    std::string msg;

    genomeErr(std::string msg_) { msg = msg_; }
};

genome_t g_empty();
genome_t g_fromfile(std::string pedname, std::string mapname);

// number of individuals
int g_nsample(genome_t g);

// number of SNPs
int g_nsnp(genome_t g);

// Removes all SNPs from genome g not present in filt
void g_filterby(genome_t g, genome_t filt);

// Filters out genomes not present among n given people in array ids
void g_filterindiv(genome_t g, std::string* ids, int n);

// Humans have 22 autosomes (present in everyone, two copies each) plus sex
// chromosomes and mitochondrial DNA. That's complex. For a first pass, remove
// everything not an autosome (chromosomes 1-22).
void g_filterchrom(genome_t g, int chromosome);

// Lookup SNP list by person
snp_t* g_plookup(genome_t g, std::string pid);

// Lookup SNP by person and SNP identifier
snp_t g_sidlookup(genome_t g, std::string pid, std::string sid);

// Lookup SNP by person and bp index
snp_t g_indlookup(genome_t g, std::string pid, int ind);

// Gets the genetic distance between SNPs i and i+1
double g_rec_dist(genome_t g, int i);

bool s_query(snp_t s, allele which);

#endif

