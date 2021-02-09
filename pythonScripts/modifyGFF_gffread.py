#! usr/bin/env/ python
import sys

GFF_in = open(sys.argv[1], 'r')
GFF_out = open(sys.argv[2], 'w')

for line in GFF_in:
    div_line = line.split("\t")
    attributes = div_line[8]
    attributes = attributes.split(";")
    new_attributes = attributes[0] + ";" + attributes[1] + "_" + attributes[5]
    div_line[8] = new_attributes
    GFF_out.write("\t".join(div_line) + "\n")

GFF_in.close()
GFF_out.close()

