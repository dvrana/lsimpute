
// DON'T TOUCH THESE -- needed to make test suite work
#include "testproto.h"

// extra libs, etc
#include <functional>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>

auto FNMAP = std::map<std::string, std::function<void(std::string)>>();

// Functions that perform the actual unit testing

void testSingleGenomeLoad(std::string args) {
    std::string pname, map;
    auto stream = std::stringstream(args);

    stream >> pname >> map;

    genome_t g = g_fromfile(pname, map);

    // TODO: actually format the output somehow
}

/* =====================================================================
 *
 * Register your functions here! The format is
 *
 *   FNMAP["tag"] = &functionname
 */
void registerFns() {
    FNMAP["loadsinglegenome"] = &testSingleGenomeLoad;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "usage: test fntag arg\n";
    }

    registerFns();
    auto fnname = std::string(argv[1]);

    FNMAP[fnname](argv[2]);
    return 0;
}

