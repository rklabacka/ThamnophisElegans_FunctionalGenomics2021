#! usr/bin/env/ python
import sys

GFF_read = open(sys.argv[1], 'r')
Variants_read = open(sys.argv[2], 'r')
TranscriptLengths = open(sys.argv[3], 'w')
log = open('log.txt', 'w')

class Gene(object):
    def __init__(self, name_in):
        self.name = name_in
        self.leng_list = []
        self.nvar = 0
        self.longestCDS = 0
    # getters
    def getName(self):
        return self.name
    def getLengs(self):
        return self.leng_list
    def getNvar(self):
        return self.nvar
    def getLongestCDS(self):
        return self.longestCDS
    # setters
    def setNvar(self, nvar_in):
        self.nvar = nvar_in
    def addLeng(self, leng_in):
        self.leng_list.append(leng_in)
    def setLongestCDS(self):
        self.longestCDS = max(self.leng_list)
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

gene_dict[past_gene].setLongestCDS()

for line in Variants_read:
    line = line.split(" ")
    haystack_gene = line[0]
    nvar = line[1].strip()
    for needle in gene_dict:
        log.write("needle: " + gene_dict[needle].getName())
        if gene_dict[needle].getName().strip() == haystack_gene.strip():
            gene_dict[needle].setNvar(nvar)
            log.write("\tNvar: " + str(gene_dict[needle].getNvar()) + ", should be: " + str(nvar))
            log.write("\tLongestCDS: " + str(gene_dict[needle].getLongestCDS()))
            propVar = float(gene_dict[needle].getNvar()) / float(gene_dict[needle].getLongestCDS())
            TranscriptLengths.write(gene_dict[needle].getName() + "\t" + str(gene_dict[needle].getNvar()) + "\t" + str(gene_dict[needle].getLongestCDS()) + "\t" + str(propVar)+ "\n")


GFF_read.close()
Variants_read.close()
TranscriptLengths.close()
