
# infrastructure
DEBUG=1

OPT=O3
CC=g++
NVCC=nvcc
CFLAGS=-std=c++11 -$(OPT)
NVCCFLAGS=-$(OPT) -m64 --gpu-architecture compute_61 -std=c++11
OBJDIR=objs
SRCDIR=src
TESTDIR=tests

LDFLAGS=-L/usr/local/depot/cuda-8.0/lib64/ -lcudart

# Not particularly important, but useful if code structure changes
PLINK=plinker
LS=hmm
EXECUTABLE=lsimpute
MAIN=$(SRCDIR)/$(EXECUTABLE).cpp

PLINKDIR=$(SRCDIR)/$(PLINK)
PLINKER=$(OBJDIR)/$(PLINK).o

HMMDIR=$(SRCDIR)/$(LS)
HMM=$(OBJDIR)/$(LS).o

LSIMPUTE_CU=lsimpute

HEADERS=$(PLINKDIR)/genome_c.h $(HMMDIR)/ls.h $(SRCDIR)/$(LSIMPUTE_CU).h

TEST_EX_NAME=tests
TEST_EX=$(TESTDIR)/$(TEST_EX_NAME)
TEST_SCRIPT=tester.py

# For every distinct "module", there should be an entry here.
OBJS=$(OBJDIR)/$(PLINK).o $(OBJDIR)/$(LS).o $(OBJDIR)/$(LSIMPUTE_CU).o

.PHONY: dirs clean debug runtest

$(EXECUTABLE): dirs $(OBJS) $(MAIN)
	$(CC) $(CFLAGS) $(LDFLAGS) -DDEBUG=0 -o $@ $(OBJS) $(MAIN)

dirs:
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(EXECUTABLE) $(OBJDIR) $(TEST_EX)

debug: DEBUG=1
debug: $(EXECUTABLE) $(TEST_EX)

# For each distinct "module", there should be a rule here. For the most part,
# the dependencies should be only the source and header files associated with
# a given module.
$(PLINKER): $(PLINKDIR)/genome.cpp $(PLINKDIR)/genome_c.h
	$(CC) $< $(CFLAGS) -c -o $@ -DDEBUG=$(DEBUG)

$(OBJDIR)/$(LSIMPUTE_CU).o: $(SRCDIR)/$(LSIMPUTE_CU).cu $(SRCDIR)/$(LSIMPUTE_CU).h
	$(NVCC) $< $(NVCCFLAGS) -c -o $@ -DDEBUG=$(DEBUG)

$(HMM): $(HMMDIR)/ls.c $(HMMDIR)/ls.h $(PLINKDIR)/genome_c.h
	$(CC) $< $(CFLAGS) -c -o $@ -DDEBUG=$(DEBUG)

$(OBJS): dirs

# Testing infrastructure

$(TEST_EX): $(OBJS)
	cd $(TESTDIR) && $(MAKE) $(TEST_EX_NAME)

runtest: debug
	cd $(TESTDIR) && ./$(TEST_EX_NAME)

