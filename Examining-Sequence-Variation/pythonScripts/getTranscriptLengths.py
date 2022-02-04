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
    argv[4] : log file
"""

import sys
import re


class Gene(object):
    def __init__(self, name_in):
        self.name = name_in
        self.leng_list = []
        self.nvar = 0
        self.pervar = 0
        self.longestCDS = 0
    def getName(self):
        return self.name
    def getLeng(self):
        return self.leng
    def getNvar(self):
        return self.nvar
    def getPervar(self):
        return self.pervar
    def getLongestCDS(self):
        return self.longestCDS
    # setters
    def setNvar(self, nvar_in):
        self.nvar = nvar_in.strip()
        if self.longestCDS != 0:
            self.pervar = float(self.nvar) / float(self.longestCDS)
    def addLeng(self, leng_in):
        self.leng_list.append(leng_in)
    def setLongestCDS(self):
        self.longestCDS = max(self.leng_list)

def get_gene_dict(gff_in):
    gene_dict = {}
    leng = 0
    past_gene = "start"
    past_cds = "start"    
    for line in GFF_read:
        line = line.split("\t")
        curr_gene = line[8]
        curr_gene = curr_gene.split("gene=")[1]
        curr_gene = curr_gene.split(";")[0]
        curr_cds = line[8].split("ID=cds-")[1]
        curr_cds = curr_cds.split(";")[0]
        if curr_cds != past_cds:
            if past_gene in gene_dict:
                gene_dict[past_gene].addLeng(leng)
            leng = 0
            past_cds = curr_cds
            if curr_gene != past_gene:
                if past_gene in gene_dict:
                    gene_dict[past_gene].setLongestCDS()
                    log.write("\nGene: " + gene_dict[past_gene].getName() + ", longest CDS: " + str(gene_dict[past_gene].getLongestCDS()))
                gene_obj = Gene(curr_gene)
                gene_dict[curr_gene] = gene_obj
                past_gene = curr_gene
        leng = leng + (int(line[4]) - int(line[3]))    
    return gene_dict

def get_variant_lengths(gene_dict, variants_read, transcript_lengths):
    for line in variants_read:
        re_string = r'total variants.+$'
        re_obj = re.compile(re_string)
        if re_obj.search(line):
            break
        line = line.split("  ")
        haystack_gene = line[0]
        nvar = line[1].strip()
        for needle in gene_dict:
            if gene_dict[needle].getName().strip() == haystack_gene.strip():
                gene_dict[needle].setNvar(nvar)
                transcript_lengths.write(gene_dict[needle].getName() + "\t" + str(gene_dict[needle].getNvar()) + "\t" + str(gene_dict[needle].getLongestCDS()) + "\t" + str(gene_dict[needle].getPervar()) + "\n")

if __name__ == "__main__":

    log = open(sys.argv[4], 'w')
    gene_dict = {}

    with open(sys.argv[1], 'r') as GFF_read:
        gene_dict = get_gene_dict(GFF_read)
        with open(sys.argv[2], 'r') as variants_read:
            with open(sys.argv[3], 'w') as transcript_lengths:
                get_variant_lengths(gene_dict, variants_read, transcript_lengths)
    log.close()

