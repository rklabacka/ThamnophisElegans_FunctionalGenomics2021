#! usr/bin/env/ python
import sys
import copy
from collections import OrderedDict

GFF_read = open(sys.argv[1], 'r')
GFF_write = open(sys.argv[2], 'w')


for line in GFF_read:
    parts = line.split("\t")
    if parts[2] == "gene":
        GFF_write.write(line)

GFF_read.close()
GFF_write.close()
