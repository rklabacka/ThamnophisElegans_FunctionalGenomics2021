#! usr/bin/env/ python
''' Script to create one file  
with all population information '''

#### Imports
import sys
import re
############

def open_file_list(file_list, outstream):
    "open file list and output content to outstream"
    with open(file_list, 'r') as instream:
        for line in instream:
            open_pop_file(line.strip(), outstream)
            
def open_pop_file(pop_file, outstream):
    "open file and copy contents to outfile"
    pop = pop_file.split(".")[0]
    with open(pop_file, 'r') as instream:
        for line in instream:
            outstream.write(line.strip() + "\t" + pop + "\n")

if __name__ == '__main__':
    files = sys.argv[1]
    output = open(sys.argv[2], 'w')
    open_file_list(files, output)
    output.close()

