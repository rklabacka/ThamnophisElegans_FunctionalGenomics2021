#! /tools/python-3.5.2/bin/ python
import sys
import re
import copy
from collections import OrderedDict
from collections import Counter

GFF_read = open(sys.argv[1], 'r')
foundGenes = open(sys.argv[2], 'r')
GFF_genes_out = open(sys.argv[3], 'w')
GFF_exons_out = open(sys.argv[4], 'w')
GFF_cds_out = open(sys.argv[5], 'w')
log = open(sys.argv[6], 'w')

foundGenes_list = []

for line in foundGenes:
    foundGenes_list.append(line.strip())

for region in GFF_read:
    region_split = region.split("\t")
    kind = region_split[2]
    if "gene=" in region_split[8]:
        gene_inquire = region_split[8].split("gene=")[1]
        gene_inquire = gene_inquire.split(";")[0]
        log.write("GFF gene: " + gene_inquire + "\n")
        if gene_inquire in foundGenes_list:
            if kind == "gene":
                GFF_genes_out.write(region) 
            elif kind == "exon":
                GFF_exons_out.write(region)
            elif kind == "CDS":
                GFF_cds_out.write(region)
        else:
            log.write(gene_inquire + " has no match\n")

GFF_read.close()
foundGenes.close()
GFF_genes_out.close()
GFF_exons_out.close()
GFF_cds_out.close()
log.close()
