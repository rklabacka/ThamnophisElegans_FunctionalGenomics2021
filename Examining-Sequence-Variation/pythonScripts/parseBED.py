#! usr/bin/env/ python
import sys

BED_in = open(sys.argv[1], 'r')
GeneList_out = open(sys.argv[2], 'w')
region = sys.argv[3]

for line in BED_in:
    line_split = line.split("\t")
    ID = line_split[3]
    ID = ID.split("-")[1]
    # Below the file name is appended with _Exons because that is what we are currently focused on. It can be changed if needed.
    outfile = ID + "_" + region  + ".bed"
    with open(outfile,'w') as BED_out:
       BED_out.write(line) 
    GeneList_out.write(ID + "\n")

BED_in.close()
GeneList_out.close()

