
#include "genome_c.h"

#include <memory>
#include <new>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <exception>

#define IFCHK(cond,onfail) if (!(cond)) { onfail; }
#define ERROR(fname,ln,msg) std::cerr << (fname) << "." << (ln) << ":" << msg

static inline int select(std::string s, std::string fname, int i) {
    switch (s[0]) {
      case 'A':
        return 0;
      case 'C':
        return 1;
      case 'G':
        return 2;
      case 'T':
        return 3;
      default:
        ERROR(fname, i, "invalid chromosome: " << s << std::endl);
        return -1;
    }
}

static inline int rselect(allele a) {
    switch(a) {
      case A:
        return 0;
      case C:
        return 1;
      case G:
        return 2;
      case T:
        return 3;
    }
}

int nsample(genome_t g) { return g->nsample; }
int nsnp(genome_t g) { return (g->map).nsnp; }

void g_filterby(genome_t g, genome_t f) {
    auto gm = g->map;
    auto fm = f->map;

    // TODO
}

void g_filterindiv(genome_t g, std::string* ids, int n) {
    auto removes = std::vector<std::string>();
    for (auto kv : g->samples) {
        if (find(ids, ids+n, kv.first) == ids+n) {
            removes.push_back(kv.first);
        }
    }

    for (auto n : removes) {
        g->nsample -= 1;
        (g->samples).erase(n);
    }
}

void g_filterchrom(genome_t g, int chromosome) {
  // TODO
}

snp_t* g_plookup(genome_t g, std::string pid) {
    auto lst = (g->samples).find(pid);

    // not found
    if (lst == (g->samples).end()) { return NULL; }

    auto result = new snp_t[(g->map).nsnp];

    for (int i = 0 ; i < (g->map).nsnp ; i += 1) {
        result[i] = lst->second[i];
    }

    return result;
}

snp_t g_sidlookup(genome_t g, std::string pid, std::string sid) {
    auto snps = g_plookup(g, pid);

    if (snps == NULL) { throw genomeErr("pid not found"); }

    auto index = (g->map).sids.find(sid);
    if (index == (g->map).sids.end()) { throw genomeErr("sid not found"); }

    return snps[index->second];
}

snp_t g_indlookup(genome_t g, std::string pid, int ind) {
    auto snps = g_plookup(g, pid);

    if (snps == NULL) { throw genomeErr("pid not found"); }

    return snps[ind];
}

bool s_query(snp_t s, allele which) {
    return s.test(rselect(which));
}

genome_t g_empty() {
    auto result = std::shared_ptr<struct genome>(new struct genome);
    (result->map).nsnp = 0;
    result->nsample = 0;
    return result;
}

// TODO: better error checking
genome_t g_fromfile(std::string pedname, std::string mapname) {
    auto result = std::shared_ptr<struct genome>(new struct genome);
    (result->map).nsnp = -1;
    result->nsample = 0;

    std::ifstream ped;
    std::ifstream map;

    std::string line;
    int ln = 0;

    map.open(mapname);
    auto data =
        std::shared_ptr<std::vector<struct snpmeta>>
        (new std::vector<struct snpmeta>);

    while (std::getline(map, line)) {
        struct snpmeta s;
        std::string chs;
        std::stringstream lstr(line);
        IFCHK(lstr >> chs,
                ERROR(mapname, ln, "parse error\n");
                return NULL);
        if (chs.compare("X") == 0) { s.num = 23; }
        else if (chs.compare("Y") == 0) { s.num = 24; }
        else {
            try {
                s.num = std::stoi(chs);
                if (s.num < 0 || s.num > 24) {
                    ERROR(mapname, ln, "chromosome number must be 1-22,X,Y");
                    return NULL;
                }
            }
            catch (std::invalid_argument) {
                ERROR(mapname, ln, "parse error\n");
                return NULL;
            }
        }

        IFCHK(lstr >> s.id,
                ERROR(mapname, ln, "parse error\n");
                return NULL);
        IFCHK(lstr >> s.dist,
                ERROR(mapname, ln, "parse error\n");
                return NULL);
        s.ind = ln;

        (*data).push_back(s);
        ln += 1;
    }

    int nsnp = ln;

    std::sort((*data).begin(), (*data).end());
    auto ids =
        std::shared_ptr<int[]>(
                new int[ln],
                std::default_delete<int[]>());

    for (int i = 0 ; i < ln ; i += 1) {
        auto s = (*data)[i];
        ids[s.ind] = i;
        (result->map).sids.insert(std::make_pair(s.id, i));
    }

    (result->map).nsnp = nsnp;
    (result->map).ids = ids;
    (result->map).data = data;

    ped.open(pedname);
    // XXX -- this loads the entire line into memory. Possibly look into
    // streaming word by word
    ln = 0;
    while (std::getline(ped, line)) {
        ln += 1;
        std::stringstream lstr(line);
        std::string fid;
        int iid, ptid, mtid, sx, ptype;
        IFCHK(lstr >> fid,
                ERROR(pedname, ln, "parse error\n");
                return NULL);
        IFCHK(lstr >> iid,
                ERROR(pedname, ln, "parse error\n");
                return NULL);
        IFCHK(lstr >> ptid,
                ERROR(pedname, ln, "parse error\n");
                return NULL);
        IFCHK(lstr >> mtid,
                ERROR(pedname, ln, "parse error\n");
                return NULL);
        IFCHK(lstr >> sx,
                ERROR(pedname, ln, "parse error\n");
                return NULL);
        IFCHK(lstr >> ptype,
                ERROR(pedname, ln, "parse error\n");
                return NULL);

        std::stringstream name;
        name << fid << "_" << iid;

        auto smp = std::shared_ptr<snp_t[]>(new snp_t[nsnp]);
        for (int i = 0 ; i < nsnp ; i += 1) {
            std::string a1, a2;

            smp[i] = snp_t(0);
            lstr >> a1 >> a2;

            int j;

            if ((j = select(a1, pedname, ln)) == -1) { return NULL; }
            smp[i].set(j);
            if ((j = select(a2, pedname, ln)) == -1) { return NULL; }
            smp[i].set(j+4);
        }

        (result->samples).insert(std::make_pair(name.str(), smp));
    }

    result->nsample = ln;

    // NOTE: fstreams don't need to be manually closed. thanks cpp destructors
    return result;
}

