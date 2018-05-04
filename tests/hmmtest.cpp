
#include <math.h>
#include <string>

#include "../src/plinker/genome_c.h"
#include "../src/hmm/ls.h"
#include "../src/lsimpute.h"
#include "infrastructure.h"
#include "lassert.h"

#define EPSILON 0.000001 // 1e-6
#define FEQ(x,y) (x > y ? ((x - y) < EPSILON) : ((y - x) < EPSILON))

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

void printLogMat(float* A, int nrow, int ncol) {
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      fprintf(stderr,"%f\t",exp(A[i * ncol + j]));
    }
    fprintf(stderr,"\n");
  }
}

void runSeqHMMBasicTest() {
    genome_t ref = g_fromfile(std::string(PED_TEST_02), std::string(MAP_TEST_02));
    genome_t sam = g_fromfile(std::string(PED_TEST_03), std::string(MAP_TEST_03));

    fprintf(stderr, "Loaded genome successfully.\n");

    /* For the below, the forward should be proportional to
     * [0.3214285714285714, 0.3214285714285714, 0.3214285714285714, 0.03571428571428571]
     * [0.04857252084793024, 0.4371526876313721, 0.4371526876313721, 0.07712210388932551]
     * [0.01606695280077814, 0.7175140083438533, 0.07972377870487259, 0.18669526015049598]
     * [0.005144151247723374, 0.8488864336559665, 0.11913289807045467, 0.02683651702585534]
     */

    float* P =  ls(sam, std::string("03_3_1"), ref, 0.1f, 1.0f);
    int nsnp = g_nsnp(sam);
    int nref = g_nsample(ref);

    for (int i = 0; i < nsnp; i++) {
      float x = rowSum(&(P[i * nref]),nref);
      if (!FEQ(x,1.0f)) {
        fprintf(stderr, "Row %d of results sums to %f\n",i, x);
        ASSERT(false, "Row does not sum to 1!");
      }
    }

    printLogMat(P, nsnp, nref);
}

void exportBasicSeqHMMTests() {
    auto basicTest = new TestCase();
    basicTest->name = (char*)"Basic Sequential HMM Functionality";
    basicTest->run = &runSeqHMMBasicTest;

    alltests.registerTest(basicTest);
}

