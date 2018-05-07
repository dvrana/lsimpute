#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "plinker/genome_c.h"
#include "hmm/ls.h"
#include "lsimpute.h"

// Look, geneticists make some long filenames, man
#define FILENAMEMAX 1024

const char* helpstring =
"Usage: lsimpute [OPTIONS] [REF] [SAMPLE]\
Uses the Li-Stephens model to impute sample genomes to a reference panel\n\n\
Arguments -t and -g are mandatory.\n\
  -g [N]        Specify garble parameter.  Must be a float > 0.0, < 1.0\n\
  -h            Print this message\n\
  -s            Run in sequential mode (much slower)\n\
  -t [N]        Specify theta.  Must be a float\n";

void printhelp() {
  printf(helpstring);
}

int main(int argc, char *argv[]) {
  float g = -1.0;
  float theta = -1.0;

  extern char* optarg;
  int opt;
  char* ref_files, * sam_files;
  bool sequential = false;

  // Read in and handle command line arguments
  while ((opt = getopt(argc, argv, "g:t:hs")) != -1) {
    switch(opt) {
      case 'h':
        printhelp();
        return 0;
        break;

      case 't':
        if (!optarg) {
          fprintf(stderr,"Must specify value for argument -t!\n");
          return 1;
        }
        theta = (float)atof(optarg);
        if (theta <= 0.0) {
          fprintf(stderr,"Theta must have positive value\n");
          return 1;
        }
        break;

      case 'g':
        if (!optarg) {
          fprintf(stderr,"Must specify value for argument -g!\n");
          return 1;
        }
        g = (float)atof(optarg);
        if (g >= 1.0 || g <= 0.0) {
          fprintf(stderr,"g must have value on range (0,1)\n");
          return 1;
        }
        break;
      case 's':
        sequential = true;
        break;
      case '?':
        break;
    }
  }

  if (argc < 3) {
    fprintf(stderr,"Must specify reference and sample files in args\n!");
    return 1;
  }

  ref_files = argv[argc - 2];
  sam_files = argv[argc - 1];

  if (g == -1.0) {
    fprintf(stderr,"Must specify garble rate with -g!\n");
    return 1;
  }

  if (theta == -1.0) {
    fprintf(stderr,"Must specify theta with -t!\n");
    return 1;
  }

  ////////////////////////////////////////
  // Actual processing code starts here //
  ////////////////////////////////////////

  // Get genome objects
  char mapname[FILENAMEMAX];
  char pedname[FILENAMEMAX];
  int reflen = strlen(ref_files);
  int samlen = strlen(sam_files);
  if (samlen > (FILENAMEMAX + 5) || reflen > (FILENAMEMAX + 5)) {
    fprintf(stderr,"Your filenames are too long, stopping to avoid error.\n");
    return 1;
  }

  strcpy(mapname, ref_files);
  strcpy(mapname + reflen, ".map");
  strcpy(pedname, ref_files);
  strcpy(pedname + reflen, ".ped");

  genome_t ref = g_fromfile(pedname, mapname);
  printf("Reading reference panel from %s and %s...\n", mapname, pedname);

  strcpy(mapname, sam_files);
  strcpy(mapname + samlen, ".map");
  strcpy(pedname, sam_files);
  strcpy(pedname + samlen, ".ped");

  genome_t sam = g_fromfile(pedname, mapname);
  printf("Reading imputed samples from %s and %s...\n", mapname, pedname);

  // Run Li-Stephens
  float* P;
  if (!sequential) {
    // XXX Note to Cam:
    // Replace the ID with the sample ID you want to LS
    // The number on the end is (at the moment) useless
    P = runThing(ref, sam, std::string("03_03_1"), 0, g, theta);
  }
  else {
    // XXX Note to Cam: same ID caveats as above
    P = ls(sam, std::string("03_03_1"), ref, g, theta);
  }

  // TODO: impute here!
  free(P); // TODO: the imputer should free this

  return 0;
}
