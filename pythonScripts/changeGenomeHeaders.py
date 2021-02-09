#! usr/bin/env/ python
import sys

GenomeRef = open(sys.argv[1], 'r')
GenomeAlt = open(sys.argv[2], 'r')
outfile = open(sys.argv[3], 'w')
headers = []

for line in GenomeRef:
    if line[0] == ">":
        headers.append(line)

num = 0

for line in GenomeAlt:
    if line[0] == ">":
        outfile.write(headers[num])
        num += 1
    else:
        outfile.write(line)

GenomeRef.close()
GenomeAlt.close()
outfile.close()

