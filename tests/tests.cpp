
#include <cstdio>

#include "lassert.h"
#include "infrastructure.h"

#include "plinktest.h"
#include "hmmtest.h"

TestFactory alltests;

// Export individual tests here!
void setup(void) {
    exportBasicPlinkerTests();
    exportBasicSeqHMMTests();
}

int main(void) {
    setup();
    alltests.runAll();
}

bool TestFactory::registerTest(TestCase* t) {
    tests.push_back(t);
    return true;
}

void TestFactory::runAll() {
    for (auto test : tests) {
        try {
            fprintf(stderr, "Running %s test...\n", test->name);
            test->run();
            fprintf(stderr, "%s test passed.\n", test->name);
        }
        catch (AssertionFailure& a) {
            fprintf(stderr, "========== TEST FAILED ===========\n");
            fprintf(stderr, "  %s: %s\n", a.loc, a.what());
        }
    }
}

