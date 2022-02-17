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
            for line in instream:
                sample = line.split("\t")[0]
                REpattern_string = r'([A-Z]{3})'
                REpattern_obj = re.compile(REpattern_string)
                search_result = re.search(REpattern_obj, line)
                if search_result.group(0) in L:
                    outstream.write(sample + "\tL\n")
                elif search_result.group(0) in M:
                    outstream.write(sample + "\tM\n")
                else:
                    print("ERROR: No ecotype assignment given for " + sample)
            
if __name__ == '__main__':
    pop_ecotypes = {"PIK" : "L" , "MAR": "L", "ELF": "L", "MAH": "M", "PAP": "M", "NAM": "M", "MER": "L", "STO": "L", "SUM": "M", "PVM": "M", "RON": "M", "CHR": "L", "ROC": "L"}
    L = ["PIK", "MAR", "ELF", "MER", "STO", "CHR", "ROC"]
    M = ["MAH", "PAP", "NAM", "SUM", "PVM", "RON"]
    readTable(sys.argv[1], sys.argv[2])
