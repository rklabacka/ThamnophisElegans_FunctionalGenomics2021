#! usr/bin/env/ python
'''
Script to add ecotype information to
pairwise population comparisons
-----------------------------------
input: infile (table with pairwise
populations in first column

output: outfile name (table formated
the same as original, with new column
appended including comparison type)
'''

#### Imports
import sys
import re
############

pop_ecotypes = {"MAR": "L", "ELF": "L", "MAH": "M", "PAP": "M", "NAM": "M", "MER": "L", "STO": "L", "SUM": "M", "PVM": "M", "RON": "M", "CHR": "L", "ROC": "L"}

def infer_ecotype_comparison(pop1, pop2):
    '''
    Infer whether populations are the
    same or different ecotypes. 
    --------------------------------
    arguments:
    pop1 and pop2: 3 letter codes that
    match those in pop_ecotypes dictionary
    ________________________________
    Returns a string of Same or Diff
    '''
    if pop_ecotypes[pop1] == pop_ecotypes[pop2]:
        return "Same"
    else:
        return "Diff"

def readTable(infile, outfile):
    '''
    Parse intable to get populations
    Add ecotype comparison info using 
    infer_ecotype_comparison funct.
    -----------------------------
    arguments:
    infile: input table (first column should 
    contain pop_pop.txt
    outfile: name for outfile
    '''
    with open(infile, 'r') as instream:
        with open(outfile, 'w') as outstream:
            outstream.write(instream.readline().strip() + "\tEcotypeComp\n")
            for line in instream:
                pops = line.split("_")
                pop1 = pops[0]
                pop2 = pops[1].split("\t")[0]
                outstream.write(line.strip() + "\t" + infer_ecotype_comparison(pop1, pop2) + "\n")
            
if __name__ == '__main__':
    readTable(sys.argv[1], sys.argv[2])
