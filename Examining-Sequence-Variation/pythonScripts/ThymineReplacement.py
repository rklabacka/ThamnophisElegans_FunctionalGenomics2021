#! usr/bin/env python
import sys

infile = open(sys.argv[1], "r")
outfile = open(sys.argv[2], "w")

for line in infile:
  if line[0] != ">":
    line = line.replace("U","T")
  outfile.write(line)
