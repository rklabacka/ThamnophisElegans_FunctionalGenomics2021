#! usr/bin/env python
import sys

depthFile = open(sys.argv[1], "r")
outfile = open(sys.argv[2], "a")

Loci = []

for line in depthFile:
    line = line.split("\t")
    gene = line[0]
    depth = line[2]
    if int(depth) > 10:
        if gene not in Loci:
            # Add locus to Loci list
            Loci.append(gene)
            outfile.write(gene + "\n")
            print(str(line))

depthFile.close()
outfile.close()


