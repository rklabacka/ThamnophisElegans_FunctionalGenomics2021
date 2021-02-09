#! usr/bin/env/ python
import sys
import copy
from collections import OrderedDict

##### CLASS LOCUS #####
class Locus(object):
    def __init__(self, exon_in, locus_name_in, start_in, stop_in):
        self.exon = exon_in
        self.locus_name = locus_name_in
        self.start = start_in
        self.stop = stop_in
### Getters
    def getExon(self):
        return self.exon
    def getName(self):
        return self.locus_name
    def getStart(self):
        return self.start
    def getStop(self):
        return self.stop

#Open infile with all genes (blast output from BlastResults.txt). Infile nature is very important:
#   Column 0: Locus name FROM INFERRED EXONS
#   Column 1: Region name FROM GENOME
#   Column 2: e value
#   Column 3: start point FROM INFERRED EXONS
#   Column 4: stop point FROM INFERRED EXONS
#   Column 5: start point FROM GENOME
#   Column 6: stop point FROM GENOME

BlastResults = open(sys.argv[1], 'r')
GFF_read = open(sys.argv[2], 'r')
GFF_write = open(sys.argv[3], 'w')
log = open(sys.argv[4], 'w')

locus_dict = {}
num = 0

for line in BlastResults:
    line = line.split("\t")
    exon = line[0]
    exon = exon + "_" + str(num)
    locus = line[1]
    start = line[5]
    stop = line[6].strip()
    new_locus = Locus(exon, locus, start, stop)
    locus_dict[exon] = new_locus
    num += 1

gff_dict = {}
for line_unsplit in GFF_read:
    line = line_unsplit.split("\t")
    haystack_locus = line[0]
    kind = line[2]
    haystack_start = line[3]
    haystack_stop = line[4]
    if kind == "gene":
        for needle in locus_dict:
            if locus_dict[needle].getName() == haystack_locus:
                if ((int(haystack_start) <= int(locus_dict[needle].getStart())) and (int(haystack_stop) >= int(locus_dict[needle].getStop()))):
                   GFF_write.write(line_unsplit) 
                   log.write("Exon <" + locus_dict[needle].getExon() + " start: " + str(locus_dict[needle].getStart()) + ", stop: " + str(locus_dict[needle].getStop()) + "> found in GFF as: <" + line_unsplit.strip() + " ; with start: " + str(haystack_start) + " and stop: " + str(haystack_stop) + ">\n")

BlastResults.close()
GFF_read.close()
GFF_write.close()
log.close()
