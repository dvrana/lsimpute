
#include <math.h>
#include <string>

#include "../src/plinker/genome_c.h"
#include "../src/hmm/ls.h"
#include "../src/lsimpute.h"
#include "infrastructure.h"
#include "lassert.h"

#define EPSILON 0.000001 // 1e-6
#define FEQ(x,y) x > y ? ((x - y) < EPSILON) : ((y - x) < EPSILON)

const char* PED_TEST_02 = "data/02.ped";
const char* MAP_TEST_02 = "data/02.map";

const char* PED_TEST_03 = "data/03.ped";
const char* MAP_TEST_03 = "data/03.map";

// Takes a log-scaled float array n and returns its non-log sum
float rowSum(float* A, int n) {
  float x = 0.0f;
  for (int i = 0; i < n; i++) x += exp(A[i]);
  return x;
}

void runSeqHMMBasicTest() {
    genome_t ref = g_fromfile(std::string(PED_TEST_02), std::string(MAP_TEST_02));
    genome_t sam = g_fromfile(std::string(PED_TEST_03), std::string(MAP_TEST_03));

    fprintf(stderr, "Loaded genome successfully.\n");

    float* P =  ls(sam, std::string("03_03_1"), ref, 0.1f, 1.0f);
    int nsnp = g_nsnp(sam);
    int nref = g_nsample(ref);

    for (int i = 0; i < g_nsnp(sam); i++) {
      float x = rowSum(&(P[i * nref]),nref);
      if (!FEQ(x,1.0f)) {
        fprintf(stderr, "Row %d of results sums to %f\n",i, x);
        ASSERT(false, "Row does not sum to 1!");
      }
    }
}

void exportBasicSeqHMMTests() {
    auto basicTest = new TestCase();
    basicTest->name = (char*)"Basic Sequential HMM Functionality";
    basicTest->run = &runSeqHMMBasicTest;

    alltests.registerTest(basicTest);
}

