#! usr/bin/env/ python
"""
Script to get: Length of transcripts
_________________________________
Arguments:
    argv[1] : gff with genes of interest
    argv[2] : tab-delimited text file with
              gene names (column 1) and
              number of variant sites
              (column 2)
    argv[3] : title of output file
"""

import sys


class Gene(object):
    def __init__(self, name_in, leng_in):
        self.name = name_in
        self.leng = leng_in
        self.nvar = 0
        self.pervar = 0
    def getName(self):
        return self.name
    def getLeng(self):
        return self.leng
    def getNvar(self):
        return self.nvar
    def getPervar(self):
        return self.pervar
    def setNvar(self, nvar_in):
        self.nvar = nvar_in.strip()
        self.pervar = float(self.nvar) / float(self.leng)

def get_gene_dict(gff_in):
    gene_dict = {}
    leng = 0
    past_gene = "start"
    
    for line in GFF_read:
        line = line.split("\t")
        curr_gene = line[8]
        curr_gene = curr_gene.split("gene=")[1]
        curr_gene = curr_gene.split(";")[0]
        if curr_gene != past_gene:
            gene_obj = Gene(past_gene, leng)
            gene_dict[past_gene] = gene_obj
            leng = 0
            past_gene = curr_gene
        leng = leng + (int(line[4]) - int(line[3]))    
    return gene_dict

def get_variant_lengths(gene_dict, variants_read):
    for line in variants_read:
        line = line.split(" ")
        haystack_gene = line[0]
        nvar = line[1]
        for needle in gene_dict:
            if gene_dict[needle].getName().strip() == haystack_gene.strip():
                gene_dict[needle].setNvar(nvar)
                transcript_lengths.write(gene_dict[needle].getName() + "\t" + str(gene_dict[needle].getNvar()) + "\t" + str(gene_dict[needle].getLeng()) + "\t" + str(gene_dict[needle].getPervar()) + "\n")

if __name__ == "__main__":

    gene_dict = {}

    with open(sys.argv[1], 'r') as GFF_read:
        gene_dict = get_gene_dict(GFF_read)
        with open(sys.argv[2], 'r') as variants_read:
            with open(sys.argv[3], 'w') as transcript_lengths:
                get_variant_lengths(gene_dict, variants_read, transcript_lengths)


