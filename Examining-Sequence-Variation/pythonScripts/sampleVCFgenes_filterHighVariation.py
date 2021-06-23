#! usr/bin/env/ python
import sys
import copy
import random
from collections import OrderedDict

##### CLASS LOCUS #####
class Gene(object):
    def __init__(self, name_in, firstSite):
        self.geneName = name_in
        self.siteList = [firstSite]
	self.length = 0
### Accessors
    def getGeneName(self):
        return self.geneName
    def getSites(self):
        return self.siteList
    def getLength(self):
        return self.length
### Mutators
    def addSite(self, nextSite):
        self.siteList.append(nextSite)
    def addLength(self, length_in):
        self.length = length_in


#Open infile with all genes (blast output from BlastResults.txt). Infile nature is very important:
#   Column 0: Gene name FROM INFERRED EXONS
#   Column 1: Region name FROM GENOME
#   Column 2: e value
#   Column 3: start point FROM INFERRED EXONS
#   Column 4: stop point FROM INFERRED EXONS
#   Column 5: start point FROM GENOME
#   Column 6: stop point FROM GENOME

VCF_read = open(sys.argv[1], 'r')
GFF = open(sys.argv[2], 'r')
VCF_write = open(sys.argv[3], 'w')
log = open(sys.argv[4], 'w')

trigger = "GENE="
geneDict = {}
count = 1 

for line in VCF_read:
    log.write("line" + str(count) + "\n")
    if "GENE=" in line:
        name = line.split("GENE=",1)[1]
        name = name.split(";")[0]
	name = name.split("-")[1]
	log.write("GENE NAME: " + name+ "\n")
        if name in geneDict:
            log.write(name + " already in geneDict\n")
            geneDict[name].addSite(line)
        else:
            log.write(name + " not in geneDict. Adding now...\n")
            new_gene = Gene(name, line)
            geneDict[name] = new_gene
    elif line.startswith("#"):
        VCF_write.write(line)
    else:
        log.write("No \"GENE=\" in line\n")
    count = count + 1

log.write("\nGeneList:\n")

for line in GFF:
    line = line.split("\t")
    gene_start = int(line[3])
    gene_end = int(line[4])
    gene_name = line[8]
    gene_length = gene_end - gene_start
    gene_name = gene_name.split("=")[1]
    gene_name = gene_name.split(";")[0]
    gene_name = gene_name.split("-")[1]
    if gene_name in geneDict:
        geneDict[gene_name].addLength(gene_length)

for gene in geneDict:
    log.write(geneDict[gene].getGeneName() + "\n")
    siteList = geneDict[gene].getSites()
    numVarSites = len(siteList)
    numTotSites = geneDict[gene].getLength()
    VarFreq = float(numVarSites) / float(numTotSites)
    log.write("\tnumber of variable sites in " + geneDict[gene].getGeneName() + ": " + str(numVarSites) + "\n")
    log.write("\tnumber of total sites in " + geneDict[gene].getGeneName() + ": " + str(numTotSites) + "\n")
    if VarFreq < 0.01:
        log.write("\tVariation frequency appropriate: " + str(VarFreq) + "\n")
        randSiteNum = random.randint(0, (numVarSites - 1))
        chosenSite = siteList[randSiteNum]
        log.write("\tchosen site: " + str(randSiteNum) + chosenSite + "\n")
        VCF_write.write(chosenSite)
    else:
        log.write("\tVariation frequency too high: " + str(VarFreq) + "\n")

VCF_read.close()
GFF.close()
VCF_write.close()
log.close()
