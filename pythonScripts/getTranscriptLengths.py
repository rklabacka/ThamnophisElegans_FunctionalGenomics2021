#! usr/bin/env/ python
import sys

GFF_read = open(sys.argv[1], 'r')
Variants_read = open(sys.argv[2], 'r')
TranscriptLengths = open(sys.argv[3], 'w')

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

for line in Variants_read:
    line = line.split(" ")
    haystack_gene = line[0]
    nvar = line[1]
    for needle in gene_dict:
        if gene_dict[needle].getName().strip() == haystack_gene.strip():
            gene_dict[needle].setNvar(nvar)
            TranscriptLengths.write(gene_dict[needle].getName() + "\t" + str(gene_dict[needle].getNvar()) + "\t" + str(gene_dict[needle].getLeng()) + "\t" + str(gene_dict[needle].getPervar()) + "\n")


GFF_read.close()
Variants_read.close()
TranscriptLengths.close()
