
# infrastructure
DEBUG=0
CC=g++
CFLAGS=-std=c++11 -DDEBUG=$(DEBUG)
OBJDIR=objs
SRCDIR=src
TESTDIR=tests

# Not particularly important, but useful if code structure changes
PLINK=plinker
LS=hmm
EXECUTABLE=lsimpute
MAIN=$(SRCDIR)/$(EXECUTABLE).cpp

PLINKDIR=$(SRCDIR)/$(PLINK)
PLINKER=$(OBJDIR)/$(PLINK).o

HMMDIR=$(SRCDIR)/$(LS)
HMM=$(OBJDIR)/$(LS).o

# Used by the testing infrastructure. Add every header file here.
HEADERS=$(PLINKDIR)/genome_c.h $(HMMDIR)/ls.h

TEST_H=$(TESTDIR)/testproto.h
TEST_EX=$(TESTDIR)/test
TEST_SCRIPT=tester.py

# For every distinct "module", there should be an entry here.
OBJS=$(OBJDIR)/$(PLINK).o $(OBJDIR)/$(LS).o

.PHONY: dirs clean runtests

$(EXECUTABLE): dirs $(OBJS) $(MAIN)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(MAIN)

dirs:
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(EXECUTABLE) $(OBJDIR) $(TEST_EX)

# For each distinct "module", there should be a rule here. For the most part,
# the dependencies should be only the source and header files associated with
# a given module.
$(PLINKER): $(PLINKDIR)/genome.cpp $(PLINKDIR)/genome_c.h
	$(CC) $< $(CFLAGS) -c -o $@

$(HMM): $(HMMDIR)/ls.c $(HMMDIR)/ls.h $(PLINKDIR)/genome_c.h
	$(CC) $< $(CFLAGS) -c -o $@

$(OBJS): dirs

# Testing infrastructure

$(TEST_EX): $(TEST_H) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(TEST_EX).cpp

$(TEST_H): $(HEADERS)
	cd $(TESTDIR) && python3 gentesth.py $(HEADERS)

runtests: $(TEST_EX)
	cd $(TESTDIR) && python3 $(TEST_SCRIPT)

