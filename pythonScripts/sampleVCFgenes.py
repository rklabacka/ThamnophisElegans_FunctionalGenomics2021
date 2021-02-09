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
### Accessors
    def getGeneName(self):
        return self.geneName
    def getSites(self):
        return self.siteList
### Mutators
    def addSite(self, nextSite):
        self.siteList.append(nextSite)


#Open infile with all genes (blast output from BlastResults.txt). Infile nature is very important:
#   Column 0: Gene name FROM INFERRED EXONS
#   Column 1: Region name FROM GENOME
#   Column 2: e value
#   Column 3: start point FROM INFERRED EXONS
#   Column 4: stop point FROM INFERRED EXONS
#   Column 5: start point FROM GENOME
#   Column 6: stop point FROM GENOME

VCF_read = open(sys.argv[1], 'r')
VCF_write = open(sys.argv[2], 'w')
log = open(sys.argv[3], 'w')

trigger = "GENE="
geneDict = {}
count = 1 

for line in VCF_read:
    log.write("line" + str(count) + "\n")
    if "GENE=" in line:
        name = line.split("GENE=",1)[1]
        name = name.split(";")[0]
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

for gene in geneDict:
    log.write(geneDict[gene].getGeneName() + "\n")
    siteList = geneDict[gene].getSites()
    numSites = len(siteList)
    log.write("\tnumber of sites in " + geneDict[gene].getGeneName() + ": " + str(numSites) + "\n")
    randSiteNum = random.randint(0, (numSites - 1))
    chosenSite = siteList[randSiteNum]
    log.write("\tchosen site: " + str(randSiteNum) + chosenSite + "\n")
    VCF_write.write(chosenSite)

VCF_read.close()
VCF_write.close()
log.close()
