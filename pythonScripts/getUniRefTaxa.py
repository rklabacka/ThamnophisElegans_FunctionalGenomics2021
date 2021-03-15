#! usr/bin/env/ python
import sys
import re

Fasta_in = open(sys.argv[1], 'r')
Taxa_out = open(sys.argv[2], 'w')

REpattern_string = r'Tax=\w+'
REpattern_obj = re.compile(REpattern_string)

for line in Fasta_in:
    if REpattern_obj.search(line):
        Taxa_out.write(REpattern_obj.search(line).group())

Fasta_in.close()
Taxa_out.close()

