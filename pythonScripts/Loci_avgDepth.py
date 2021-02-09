#! usr/bin/env python
import sys

count = 0
new_locus= ""
old_locus = ""
curr_bp = 0
old_bp = 0
depthSum = 0
depthAvg = 0
bpCount = 0

depthFile = open(sys.argv[1], "r")
outfile = open(sys.argv[2], "a")

for line in depthFile:
    line = line.split("\t")
    new_locus = line[0]
    curr_bp = int(line[1])
    if count > 0:
        if new_locus != old_locus:
            # new locus, calculate and report previous exon
            depthAvg = float(depthSum / bpCount)
            outfile.write(str(bpCount) + "\t" + str(depthAvg) + "\n")
            bpCount = 0
            depthSum = int(line[2])
        else:
            if (curr_bp - old_bp) > 300:
                # new exon, calculate and report previous exon
                depthAvg = float(depthSum / bpCount)
                outfile.write(str(bpCount) + "\t" + str(depthAvg) + "\n")
                bpCount = 0
                depthSum = int(line[2])
            else:
                # same exon, add to sum
                depthSum = depthSum + int(line[2])
    old_bp = curr_bp
    bpCount = bpCount + 1
    old_locus = new_locus
    count = count + 1

depthFile.close()
outfile.close()


