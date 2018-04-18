
DEBUG=0
CC=g++
CFLAGS=-std=c++11 -DDEBUG=$(DEBUG)
OBJDIR=objs
SRCDIR=src

PLINK=plinker
EXECUTABLE=lsimpute

OBJS=$(OBJDIR)/$(PLINK).o

.PHONY: dirs clean

$(EXECUTABLE): dirs $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS)

dirs:
	mkdir -p $(OBJDIR)

$(OBJDIR)/$(PLINK).o: $(SRCDIR)/$(PLINK)/genome.cpp $(SRCDIR)/$(PLINK)/genome_c.h
	$(CC) $< $(CFLAGS) -c -o $@

clean:
	rm -rf $(EXECUTABLE) $(OBJDIR)

