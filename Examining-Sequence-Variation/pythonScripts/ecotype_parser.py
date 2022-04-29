#! usr/bin/env/ python
''' Script to parse individuals into
files for respective populations '''

#### Imports
import sys
import re
############

def sort_by_ecotype(instream, re_L_obj, re_M_obj):
    '''
    Loop through instream and sort each sample by ecotype
    Then, create text files with each sample for the two ecotypes
    '''
    outfile_L = ""
    outfile_M = ""
    for line in in_stream:
        if re_L_obj.search(line):
            outfile_L += line 
        elif re_M_obj.search(line):
            outfile_M += line 
        else:
            print("ERROR: " + line.strip() + " not found")
    return outfile_L, outfile_M 
                

if __name__ == '__main__':
    pop_set = set()
    # Establish regular expression pattern for populations
    re_L_string = r'(?:MAR)|(?:ELF)|(?:PIK)|(?:MER)|(?:STO)|(?:CHR)|(?:ROC)'
    re_M_string = r'(?:MAH)|(?:PAP)|(?:NAM)|(?:SUM)|(?:PVM)|(?:RON)'
    # Compile regular expression pattern into regex object
    re_L_obj = re.compile(re_L_string)
    re_M_obj = re.compile(re_M_string)
    # Create text file for each population
    with open(sys.argv[1], 'r') as in_stream:
        outfile_L, outfile_M = sort_by_ecotype(in_stream, re_L_obj, re_M_obj)
    with open(sys.argv[2], 'w') as out_stream:
        out_stream.write(outfile_L)
    with open(sys.argv[3], 'w') as out_stream:
        out_stream.write(outfile_M)
