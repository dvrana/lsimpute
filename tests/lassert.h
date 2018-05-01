
#ifndef LASSERT_H
#define LASSERT_H
#define LASSERT_BUF_SIZE 256
#include <exception>
#include <cstdio>

struct AssertionFailure : public std::exception {
    const char* msg;
    char* loc;

    AssertionFailure(const char* msg_, int ln, const char* file) {
        msg = msg_;
        loc = new char[LASSERT_BUF_SIZE];
        std::sprintf(loc, "%s:%d", file, ln);
    }

    const char* what() const throw() {
        return msg;
    }
};

inline void _MkAssert(int x, const char* msg, int lineno, const char* file) {
    if (!x) { throw AssertionFailure(msg, lineno, file); }
}

#define ASSERT(p, msg) _MkAssert(p, msg, __LINE__, __FILE__)

#endif

