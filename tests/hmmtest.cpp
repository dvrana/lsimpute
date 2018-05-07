
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
     * [0.048532109899404895, 0.436788989094644, 0.436788989094644, 0.07788991191130712]
     * [0.01668948789147723, 0.7115906893757763, 0.07906563215286405, 0.1926541905798823]
     * [0.005474670956529986, 0.8458619051087488, 0.12077602266265874, 0.027887401272062663]
     *
     * and backward to
     * [0.03486421210153343, 0.7867259437664135, 0.15709310523136097, 0.021316738900692092]
     * [0.006624653674670816, 0.7165021824634252, 0.12202824578123218, 0.15484491808067163]
     * [0.013806503278561623, 0.7757414704929454, 0.08619349672143839, 0.1242585295070546]
     * [0.05, 0.45, 0.45, 0.05]
     *
     * P should equal
     * [0.007682004169127661, 0.8308619624778509, 0.14150498107198667, 0.019951052281034633]
     * [0.001732168759242964, 0.8759231805429016, 0.0973247978381002, 0.025019852859755082]
     * [0.0022783504437034583, 0.874279472249249, 0.09714216358324991, 0.026300013723797623]
     * [0.0054746709565299855, 0.8458619051087486, 0.12077602266265872, 0.027887401272062656]
     */

    float* P =  ls(sam, std::string("03_03_1"), ref, 0.1f, 1.0f);
    int nsnp = g_nsnp(sam);
    int nref = g_nsample(ref);

    // Test row sums
    for (int i = 0; i < nsnp; i++) {
      float x = rowSum(&(P[i * nref]),nref);
      if (!FEQ(x,1.0f)) {
        fprintf(stderr, "Row %d of results sums to %f\n",i, x);
        ASSERT(false, "Row does not sum to 1!");
      }
    }

    // Test smoothed values
    ASSERT(FEQ(exp(P[0]),0.007682004169127661),
        "Seq HMM result at (0,0) incorrect!");
    ASSERT(FEQ(exp(P[1]),0.8308619624778509),
        "Seq HMM result at (1,1) incorrect!");
    ASSERT(FEQ(exp(P[2]),0.14150498107198667),
        "Seq HMM result at (1,2) incorrect!");
    ASSERT(FEQ(exp(P[3]),0.019951052281034633),
        "Seq HMM result at (1,3) incorrect!");

    ASSERT(FEQ(exp(P[4]),0.001732168759242964),
        "Seq HMM result at (2,0) incorrect!");
    ASSERT(FEQ(exp(P[5]),0.8759231805429016),
        "Seq HMM result at (2,1) incorrect!");
    ASSERT(FEQ(exp(P[6]),0.0973247978381002),
        "Seq HMM result at (2,2) incorrect!");
    ASSERT(FEQ(exp(P[7]),0.025019852859755082),
        "Seq HMM result at (2,3) incorrect!");

    ASSERT(FEQ(exp(P[8]),0.0022783504437034583),
        "Seq HMM result at (3,0) incorrect!");
    ASSERT(FEQ(exp(P[9]),0.874279472249249),
        "Seq HMM result at (3,1) incorrect!");
    ASSERT(FEQ(exp(P[10]),0.09714216358324991),
        "Seq HMM result at (3,2) incorrect!");
    ASSERT(FEQ(exp(P[11]),0.026300013723797623),
        "Seq HMM result at (3,3) incorrect!");

    ASSERT(FEQ(exp(P[12]),0.005474670956529986),
        "Seq HMM result at (3,0) incorrect!");
    ASSERT(FEQ(exp(P[13]),0.8458619051087488),
        "Seq HMM result at (3,1) incorrect!");
    ASSERT(FEQ(exp(P[14]),0.12077602266265874),
        "Seq HMM result at (3,2) incorrect!");
    ASSERT(FEQ(exp(P[15]),0.027887401272062663),
        "Seq HMM result at (3,3) incorrect!");
}

void runGPUHMMBasicTest() {
  genome_t ref = g_fromfile(std::string(PED_TEST_02), std::string(MAP_TEST_02));
  genome_t sam = g_fromfile(std::string(PED_TEST_03), std::string(MAP_TEST_03));

  int nsample = g_nsample(ref);
  int nsnp = g_nsnp(ref);

  float* P = runThing(ref, sam, std::string("03_03_1"), 12, 0.1, 1.0);

  // Test row sums
  for (int i = 0; i < nsnp; i++) {
    float x = rowSum(&(P[i * nsample]),nsample);
    if (!FEQ(x,1.0f)) {
      fprintf(stderr, "Row %d of results sums to %f\n",i, x);
      ASSERT(false, "Row does not sum to 1!");
    }
  }

  // Test smoothed values
  ASSERT(FEQ(exp(P[0]),0.007682004169127661),
      "GPU HMM result at (0,0) incorrect!");
  ASSERT(FEQ(exp(P[1]),0.8308619624778509),
      "GPU HMM result at (1,1) incorrect!");
  ASSERT(FEQ(exp(P[2]),0.14150498107198667),
      "GPU HMM result at (1,2) incorrect!");
  ASSERT(FEQ(exp(P[3]),0.019951052281034633),
      "GPU HMM result at (1,3) incorrect!");

  ASSERT(FEQ(exp(P[4]),0.001732168759242964),
      "GPU HMM result at (2,0) incorrect!");
  ASSERT(FEQ(exp(P[5]),0.8759231805429016),
      "GPU HMM result at (2,1) incorrect!");
  ASSERT(FEQ(exp(P[6]),0.0973247978381002),
      "GPU HMM result at (2,2) incorrect!");
  ASSERT(FEQ(exp(P[7]),0.025019852859755082),
      "GPU HMM result at (2,3) incorrect!");

  ASSERT(FEQ(exp(P[8]),0.0022783504437034583),
      "GPU HMM result at (3,0) incorrect!");
  ASSERT(FEQ(exp(P[9]),0.874279472249249),
      "GPU HMM result at (3,1) incorrect!");
  ASSERT(FEQ(exp(P[10]),0.09714216358324991),
      "GPU HMM result at (3,2) incorrect!");
  ASSERT(FEQ(exp(P[11]),0.026300013723797623),
      "GPU HMM result at (3,3) incorrect!");

  ASSERT(FEQ(exp(P[12]),0.005474670956529986),
      "GPU HMM result at (3,0) incorrect!");
  ASSERT(FEQ(exp(P[13]),0.8458619051087488),
      "GPU HMM result at (3,1) incorrect!");
  ASSERT(FEQ(exp(P[14]),0.12077602266265874),
      "GPU HMM result at (3,2) incorrect!");
  ASSERT(FEQ(exp(P[15]),0.027887401272062663),
      "GPU HMM result at (3,3) incorrect!");
}

void exportBasicHMMTests() {
    auto basicTest = new TestCase();
    basicTest->name = (char*)"Basic Sequential HMM Functionality";
    basicTest->run = &runSeqHMMBasicTest;

    auto gpuTest = new TestCase();
    gpuTest->name = (char*)"Basic GPU HMM Functionality";
    gpuTest->run = &runGPUHMMBasicTest;

    alltests.registerTest(basicTest);
    alltests.registerTest(gpuTest);
}

