#! usr/bin/env/ python
"""Script used for modifying a gff file
to be compatible with gffread
"""

import sys
import re

GFF_in = open(sys.argv[1], 'r')
GFF_out = open(sys.argv[2], 'w')

for line in GFF_in:
    div_line = line.split("\t")
    attributes = div_line[8]
    pattern_str = r'(ID=[\w\-\.]+);(Parent=[\w\-\.]+);.+(gene=\w+)'
    pattern_obj = re.compile(pattern_str)
    search_result = re.search(pattern_obj, attributes)
    new_attributes = search_result.group(1) + ";" + search_result.group(2) + "_" + search_result.group(3)
    div_line[8] = new_attributes
    GFF_out.write("\t".join(div_line) + "\n")

GFF_in.close()
GFF_out.close()

