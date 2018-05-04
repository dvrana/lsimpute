
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

#if DEBUG
#include <cstdio>
#define D_PRINTF(args...) fprintf(stderr, args)
#else
#define D_PRINTF(args...)
#endif

#define IFCHK(cond,onfail) if (!(cond)) { onfail; }
#define ERROR(fname,ln,msg) std::cerr << (fname) << "." << (ln) << ":" << msg

static inline allele select(std::string s, std::string fname, int i) {
    switch (s[0]) {
      case 'A':
        return A;
      case 'C':
        return C;
      case 'G':
        return G;
      case 'T':
        return T;
      default:
        ERROR(fname, i, "invalid chromosome: " << s << std::endl);
        throw genomeErr("Unknown allele string.");
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

static inline int search(
    std::shared_ptr<std::string> arr, int n, std::string arg
) {
    std::string* arr_ = arr.get();
    for (int i = 0 ; i < n ; i += 1) {
        if (arr_[i] == arg) { return i; }
    }
    return -1;
}

int g_nsample(genome_t g) { return g->nsample; }
int g_nsnp(genome_t g) { return (g->map).nsnp; }

void g_filterby(genome_t g, genome_t f) {
    auto gm = g->map;
    auto fm = f->map;

    auto keeps = std::vector<struct snpmeta>(16);
    for (int i = 0 ; i < gm.nsnp ; i += 1) {
        // XXX - this would just be gm.id[i] if we had a more recent version of
        // g++ on andrew...
        auto idx = gm.id_arr()[i];
        auto k = (*gm.data)[idx];
        if (search(fm.sids, fm.nsnp, k.id) != -1) {
            keeps.push_back((*gm.data)[idx]);
        }
    }

    size_t nsnp = keeps.size();

    for (size_t i = 0 ; i < nsnp ; i += 1) { keeps[i].ind = i; }

    (g->map).nsnp = nsnp;
    std::sort(keeps.begin(), keeps.end());

    (g->map).sids = std::shared_ptr<std::string>(
          new std::string[nsnp],
          std::default_delete<std::string[]>()
        );
    (g->map).ids = std::shared_ptr<int>
        (new int[nsnp], std::default_delete<int[]>());

    for (size_t i = 0 ; i < nsnp ; i += 1) {
        (g->map).sids.get()[i] = keeps[i].id;
        (g->map).id_arr()[keeps[i].ind] = i;
    }

    (g->map).data = std::make_shared<std::vector<struct snpmeta>>(keeps);
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
    auto keeps = std::vector<struct snpmeta>();
    auto gm = g->map;

    for (int i = 0 ; i < gm.nsnp ; i += 1) {
        auto idx = gm.id_arr()[i];
        auto k = (*gm.data)[idx];
        if (k.chnum == chromosome) { keeps.push_back(k); }
    }

    auto nsnp = keeps.size();

    for (size_t i = 0 ; i < nsnp ; i += 1) { keeps[i].ind = i; }

    (g->map).nsnp = nsnp;
    std::sort(keeps.begin(), keeps.end());

    (g->map).ids = std::shared_ptr<int>
        (new int[nsnp], std::default_delete<int[]>());

    for (size_t i = 0 ; i < nsnp ; i += 1) {
        (g->map).sids.get()[i] = keeps[i].id;
        (g->map).id_arr()[keeps[i].ind] = i;
    }
}

double g_rec_dist(genome_t g, int index) {
    auto m = g->map.data;
    return ((*m)[index+1].gdist - (*m)[index].gdist);
}

snp_t* g_plookup(genome_t g, std::string pid) {
    auto lst = (g->samples).find(pid);

    // not found
    if (lst == (g->samples).end()) { return NULL; }

    auto result = new snp_t[(g->map).nsnp];

    for (int i = 0 ; i < (g->map).nsnp ; i += 1) {
        result[i] = (lst->second).get()[i];
    }

    return result;
}

snp_t g_sidlookup(genome_t g, std::string pid, std::string sid) {
    auto snps = g_plookup(g, pid);

    if (snps == NULL) { throw genomeErr("pid not found"); }

    auto index = search((g->map).sids, (g->map).nsnp, sid);
    if (index == -1) { throw genomeErr("sid not found"); }

    return snps[index];
}

snp_t g_indlookup(genome_t g, std::string pid, int ind) {
    auto snps = g_plookup(g, pid);

    if (snps == NULL) { throw genomeErr("pid not found"); }

    return snps[ind];
}

bool s_query(snp_t s, allele which) {
    return s == which;
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
        if (chs.compare("X") == 0) { s.chnum = 23; }
        else if (chs.compare("Y") == 0) { s.chnum = 24; }
        else {
            try {
                s.chnum = std::stoi(chs);
                if (s.chnum < 0 || s.chnum > 24) {
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
        IFCHK(lstr >> s.gdist,
                ERROR(mapname, ln, "parse error\n");
                return NULL);
        IFCHK(lstr >> s.pos,
                ERROR(mapname, ln, "parse error\n");
                return NULL);
        s.ind = ln;

        (*data).push_back(s);
        ln += 1;
    }

    int nsnp = ln;

    std::sort((*data).begin(), (*data).end());
    auto ids =
        std::shared_ptr<int>(
                new int[ln],
                std::default_delete<int[]>());
    (result->map).sids =
        std::shared_ptr<std::string>(
            new std::string[ln], std::default_delete<std::string[]>());

    for (int i = 0 ; i < ln ; i += 1) {
        auto s = (*data)[i];
        ids.get()[s.ind] = i;
        (result->map).sids.get()[i] = s.id;
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

        auto smp1 = std::shared_ptr<snp_t>(new snp_t[nsnp]);
        auto smpp1 = smp1.get();
        auto smp2 = std::shared_ptr<snp_t>(new snp_t[nsnp]);
        auto smpp2 = smp2.get();
        for (int i = 0 ; i < nsnp ; i += 1) {
            std::string a1, a2;

            lstr >> a1 >> a2;

            allele j;

            if ((j = select(a1, pedname, ln)) == -1) { return NULL; }
            smpp1[i] = j;
            if ((j = select(a2, pedname, ln)) == -1) { return NULL; }
            smpp2[i] = j;
        }

        std::stringstream n1, n2;
        n1 << name.str() << "_1";
        n2 << name.str() << "_2";
        (result->samples).insert(std::make_pair(n1.str(), smp1));
        D_PRINTF("inserting sample name %s\n", n1.str().c_str());
        (result->samples).insert(std::make_pair(n2.str(), smp2));
        D_PRINTF("inserting sample name %s\n", n2.str().c_str());
    }

    result->nsample = ln*2;

    // NOTE: fstreams don't need to be manually closed. thanks cpp destructors
    return result;
}

