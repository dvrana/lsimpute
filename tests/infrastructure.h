
#ifndef TEST_INFRASTRUCTURE
#define TEST_INFRASTRUCTURE

#include <functional>
#include <vector>

#include "lassert.h"

struct TestCase {
    char* name;
    void (*run)();
};

class TestFactory {
private:
    std::vector<TestCase*> tests;

public:
    bool registerTest(TestCase* t);
    void runAll();
};

extern TestFactory alltests;

#endif

