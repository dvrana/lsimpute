#!/usr/bin/python

import sys

PROJECT_ROOT = ".."

with open("testproto.h", 'w') as f:
    for hfile in sys.argv[1:]:
        print(f'#include "{PROJECT_ROOT}/{hfile}"', file=f)

