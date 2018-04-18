
# infrastructure
DEBUG=0
CC=g++
CFLAGS=-std=c++11 -DDEBUG=$(DEBUG)
OBJDIR=objs
SRCDIR=src

# Not particularly important, but useful if code structure changes
PLINK=plinker
EXECUTABLE=lsimpute
MAIN=$(SRCDIR)/$(EXECUTABLE).cpp

PLINKDIR=$(SRCDIR)/$(PLINK)
PLINKER=$(OBJDIR)/$(PLINK).o

# For every distinct "module", there should be an entry here.
OBJS=$(OBJDIR)/$(PLINK).o

.PHONY: dirs clean

$(EXECUTABLE): dirs $(OBJS) $(MAIN)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(MAIN)

dirs:
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(EXECUTABLE) $(OBJDIR)

# For each distinct "module", there should be a rule here. For the most part,
# the dependencies should be only the source and header files associated with
# a given module.
$(PLINKER): $(PLINKDIR)/genome.cpp $(PLINKDIR)/genome_c.h
	$(CC) $< $(CFLAGS) -c -o $@

