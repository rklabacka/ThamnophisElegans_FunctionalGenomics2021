#! usr/bin/env/ python
import sys
import copy
from collections import OrderedDict


AllExons_in = open(sys.argv[1], 'r')
WantedGenes_in = open(sys.argv[2], 'r')
WantedExons_out = open(sys.argv[3], 'w')
log = open(sys.argv[4], 'w')
WantedGenes = []
FoundGenes = []
FoundExonsCount = 0

for line in WantedGenes_in:
  WantedGenes.append(line.strip())

WantedGenes_copy = WantedGenes.copy()

log.write("WantedGenes: " + str(len(WantedGenes)) + "\n")

for line in AllExons_in:
  if line[0] == ">":
    exon=line.replace(">","")
    if exon.split("_")[0] in WantedGenes:
      if exon.split("_")[0] in WantedGenes_copy:
        WantedGenes_copy.remove(exon.split("_")[0])
      if exon.split("_")[0] not in FoundGenes:
        FoundGenes.append(exon.split("_")[0])
      FoundExonsCount += 1
      WantedExons_out.write(">" + exon)
      WantedExons_out.write(AllExons_in.readline())

log.write("FoundGenes: " + str(len(FoundGenes)) + "\n")
log.write("FoundExons: " + str(FoundExonsCount) + "\n")
log.write("Genes not found: " + str(WantedGenes_copy))

AllExons_in.close()
WantedGenes_in.close()
WantedExons_out.close()
log.close()
