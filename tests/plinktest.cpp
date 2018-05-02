
#include <string>

#include "../src/plinker/genome_c.h"
#include "infrastructure.h"
#include "lassert.h"

const char* PED_TEST_01 = "data/01.ped";
const char* MAP_TEST_01 = "data/01.map";

void runPlinkBasicTest() {
    genome_t g = g_fromfile(std::string(PED_TEST_01), std::string(MAP_TEST_01));

    fprintf(stderr, "Loaded genome successfully.\n");
    ASSERT(g_nsample(g) == 2, "genome_t should report the correct nsample");
    ASSERT(g_nsnp(g) == 3, "genome_t should report the correct nsnp");

    auto s = g_plookup(g, (char *)"01_1");
    fprintf(stderr, "Looked up sample successfully.\n");

    ASSERT(s_query(s[0], A), "first allele not correctly recorded");
    ASSERT(s_query(s[0], G), "second allele not correctly recorded");
    ASSERT(s_query(s[1], G), "second pair of alleles not recorded");
    ASSERT(s_query(s[2], A), "last pair of alleles not recorded");
    ASSERT(s_query(s[2], C), "last pair of alleles not recorded");

    auto s2 = g_plookup(g, (char*)"01_2");
    fprintf(stderr, "Looked up sample successfully.\n");

    ASSERT(s_query(s2[0], A), "failure in reading second sample");
}

void exportBasicPlinkerTests() {
    auto basicTest = new TestCase();
    basicTest->name = (char*)"Basic Plinker Functionality";
    basicTest->run = &runPlinkBasicTest;

    alltests.registerTest(basicTest);
}

