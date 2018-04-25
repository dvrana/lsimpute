
#include <vector>
#include <unordered_map>
#include <memory>
#include <string>
#include <bitset>

// SNPs are stored as bitvectors in ACGT order
typedef std::bitset<8> snp_t;

enum allele { A, C, G, T };
typedef std::pair<allele, allele> apair;

struct snpmeta {
    int ind; // index in mapfile ordering
    std::string id;
    int chnum;
    int dist;

    bool operator < (const snpmeta& s) const {
        return dist < s.dist;
    }
};

struct snpmap {
    int nsnp;
    // maps index in mapfile to index in bp ordering
    std::shared_ptr<int[]> ids;
    // maps snp id string to index in bp ordering
    std::unordered_map<std::string, int> sids;
    // ordered in bp order
    std::shared_ptr<std::vector<struct snpmeta>> data;
};

struct genome {
    int nsample;
    struct snpmap map;
    // maps familyid_indid to sample
    std::unordered_map<std::string, std::shared_ptr<snp_t[]>> samples;
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

bool s_query(snp_t s, allele which);

