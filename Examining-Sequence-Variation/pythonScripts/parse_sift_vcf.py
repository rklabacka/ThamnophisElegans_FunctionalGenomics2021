#! usr/bin/env/ python
''' Script to parse individuals into
files for respective populations '''

#### Imports
import sys
import re
############

def parse_each_line(infile_in, outfile_in):
    '''
    function that loops through infile
    and parses each line, printing the
    outfile in tab-delimitted format.
    '''
    outfile.write('CHROM\tPOS\tREF\tALT\tGENE\tAA\tAAPOS\tSIFTSCORE\tSIFTMED\tSIFTCAT\n')
    for line in infile_in:
        if line[0] != "#":
            line = line.split("\t")
            chrom = line[0]
            site = line[1]
            ref = line[3]
            alt = line[4]
            info = line[7]
            info_dict = parse_info(info)
            info = '\t'.join(info_dict.values())
            outfile.write('{}\t{}\t{}\t{}\t{}\n'.format(chrom, site, ref, alt, info))

def parse_info(info_in):
    parsed_info = {}
    pattern_str = r'\|gene-(\w+)\|.+(?:NONCODING|NONSYNONYMOUS)\|(\w\/\w)\|(\d+)\|([\d\.]+)\|([\d\.]+)\|(\d+)\|\w+\|(\w+)'
    pattern_obj = re.compile(pattern_str)
    search_result = re.search(pattern_obj, info_in)
    parsed_info["gene"] = search_result.group(1)
    parsed_info["amino_acid"] = search_result.group(2)
    parsed_info["aa_pos"] = search_result.group(3)
    parsed_info["sift_score"] = search_result.group(4)
    parsed_info["sift_median"] = search_result.group(5)
    parsed_info["n_seq"] = search_result.group(6)
    parsed_info["category"] = search_result.group(7)
    return parsed_info


if __name__ == '__main__':
    with open(sys.argv[1], 'r') as infile:
        with open(sys.argv[2], 'w') as outfile:
            parse_each_line(infile, outfile)
