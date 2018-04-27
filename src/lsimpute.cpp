
#include "lsimpute.h"
#include "plinker/genome_c.h"

#include <memory>

// TODO: clean this up
void runThing(genome_t g, std::string id) {
    auto nsnp = g_nsnp(g);
    auto thing = std::shared_ptr<lsimputer>(new lsimputer(nsnp));

    auto snpmap = g_plookup(g, id);

    for (int i = 0 ; i < nsnp ; i += 1) {
        thing->snps[i] = snpmap[i].to_ulong();
        if (i+1 >= nsnp) { thing->dists[i] = g_rec_dist(g, i); }
    }

    thing->compute();
}

int main() {
    return 0;
}

