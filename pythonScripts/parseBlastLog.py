#! usr/bin/env/ python
import sys
import copy
from collections import OrderedDict


BlastLog_in = open(sys.argv[1], 'r')
GeneBook_out = open(sys.argv[2], 'w')

GeneBook_out.write("ExonProbe	GeneName	Dbxref\n")

for line in BlastLog_in:
  Exon = line.split("<")[1]
  Exon = Exon.split(" ")[0]
  Match = line.split("Dbxref=GeneID:")[1]
  Match = Match.split(",")[0]
  Gene = line.split("gene=")[1]
  Gene = Gene.split(";")[0]
  if len(Match) > 15:
    Match = "No Annotation"
  GeneBook_out.write(Exon + "\t" + Gene + "\t" + Match + "\n")

BlastLog_in.close()
GeneBook_out.close()
