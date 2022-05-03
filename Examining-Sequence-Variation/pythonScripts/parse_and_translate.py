#! usr/bin/env/ python
import sys
from Bio.Seq import Seq

# fna_in = individual fasta file with all genes
# fna_out = fasta files for each gene of the individual
# faa_out = fasta files for each protein of the individual
fna_in = open(sys.argv[1], 'r')
IndName = sys.argv[2]
header = fna_in.readline()
fna_out = open(IndName + "_" + header.lstrip(">").rstrip("\n") + ".fna", 'a')
faa_out = open(IndName + "_" + header.lstrip(">").rstrip("\n") + ".faa", 'a')
cds = ""
pepseq = ""

for line in fna_in:
    if line[0] == ">":
        fna_out.write(">" + IndName + "_" + header.lstrip(">") + cds) 
        fna_out.close()
        Seq_cds = Seq(cds)
        pepseq = Seq_cds.translate()
        faa_out.write(">" + IndName + "_" + header.lstrip(">") + str(pepseq))
        faa_out.close()
        # Begin new isoform/gene
        header = line
        fna_out = open(IndName + "_" + header.lstrip(">").rstrip("\n") + ".fna", 'a')
        faa_out = open(IndName + "_" + header.lstrip(">").rstrip("\n") + ".faa", 'a')
        cds = ""
    else:
        cds = cds + line.strip().replace(".","N")

fna_in.close()
fna_out.close()

