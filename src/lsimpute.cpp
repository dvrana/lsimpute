
#include <memory>
#include "lsimpute.h"

#include "plinker/genome_c.h"

lsimputer::lsimputer(genome_t G, float g_, float theta_) {
    nsnp = g_nsnp(G);
    nsample = g_nsample(G);

    g = g_;
    theta = theta_;

    dists = new float[nsnp];
    ref = new uint8_t[nsnp * nsample];

    for (int i = 0 ; i < nsnp-1 ; i += 1) {
        dists[i] = g_rec_dist(G, i);
    }

    auto offs = 0;
    for (auto entry : *G) {
        auto snpmap = entry.second;
        for (int i = 0 ; i < nsnp ; i += 1) {
            ref[offs + (i * nsample)] = snpmap.get()[i];
        }
        offs += 1;
    }
}

lsimputer::~lsimputer() {
    delete ref;
    delete dists;
}

// TODO: clean this up
// TODO: make this return multiple things or reuse data structures
float* ls_gpu(genome_t G, genome_t impute, std::string id, int chr,
    float g, float theta) {
    auto thing = std::shared_ptr<lsimputer>(new lsimputer(G, g, theta));

    auto nsnp = g_nsnp(G);
    auto snpmap = g_plookup(impute, id);
    auto param = new uint8_t[nsnp];

    for (int i = 0 ; i < nsnp ; i += 1) {
        param[i] = snpmap[i];
    }

    float* P = thing->compute(param);

    delete param;
    return P;
}

