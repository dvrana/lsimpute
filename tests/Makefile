
# recursive make considered harmful xd
#
# This makefile should only be invoked via the top-level makefile.

CC=g++
OBJDIR=../objs
TOBJDIR=scratch
DEBUG=1
CFLAGS=-std=c++11 -DDEBUG=1
LDFLAGS=-L/usr/local/depot/cuda-8.0/lib64/ -lcudart

TEST_EX=tests
OBJS=$(OBJDIR)/*.o
TOBJS=$(TOBJDIR)/plinktest.o $(TOBJDIR)/hmmtest.o

.PHONY: all dirs

all:
	@echo "Don't execute make from this dir! This doesn't do proper dependency tracking!"
	@echo "Use ``make $(TEST_EX)`` from the project top level instead!"

dirs:
	mkdir -p $(TOBJDIR)

$(TEST_EX): $(TEST_EX).cpp infrastructure.h $(TOBJS)
	$(CC) $< $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(TOBJS)

$(TOBJDIR)/%.o: %.cpp %.h dirs
	$(CC) $< $(CFLAGS) -c -o $@

