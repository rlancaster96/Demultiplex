#!/usr/bin/env python

import re

f = "indexes.txt"

indexfiledict: dict = {}

with open(f, "r") as fh:
    fh.readline()
    for line in fh:
        line = line.strip()
        line = line.split()
        indexfiledict[line[4]] = line

print(indexfiledict)