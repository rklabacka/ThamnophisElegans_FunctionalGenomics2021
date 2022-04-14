#! usr/bin/env/ python
'''
Script for converting VCF with SIFT scores
into a tab-delimmited text file. An additional
VCF with SNPEFF scores can also have categories
added, but this is not required.

Parameters:
sift_annotated.vcf
(optional snpeff_annotated.vcf)
output.txt

Examples:
python parse_sift_vcf.py sift_annotated.vcf snpeff_annotated.vcf output.txt
python parse_sift_vcf.py sift_annotated.vcf output.txt
'''

#### Imports
import sys
import re
############

class AnnotatedSNP:
    '''
    Class for annotated SNPs
    '''
    
    def __init__(self, chrom_in, pos_in, ref_in, alt_in, info_in):
        self.chrom = chrom_in 
        self.pos = pos_in
        self.ref = ref_in
        self.alt = alt_in
        self.info = info_in
        self.snpeff_cat = "NA"
        self.parse_sift_info(self.info)
    
    def parse_sift_info(self, info_in):
        pattern_str = r'\|gene-(\w+)\|.+(?:NONCODING|NONSYNONYMOUS)\|(\w\/\w)\|(\d+)\|([\d\.]+)\|([\d\.]+)\|(\d+)\|\w+\|(\w+)'
        pattern_obj = re.compile(pattern_str)
        search_result = re.search(pattern_obj, info_in)
        self.gene = search_result.group(1)
        self.amino_acid = search_result.group(2)
        self.aa_pos = search_result.group(3)
        self.sift_score = search_result.group(4)
        self.sift_median = search_result.group(5)
        self.n_seq = search_result.group(6)
        self.sift_cat = search_result.group(7)

    def parse_snpeff_info(self, info_in):
        pattern_str = r'\|(?:missense_variant|missense_variant\&splice_region_variant)\|(\w+)\|'
        pattern_obj = re.compile(pattern_str)
        search_result = re.search(pattern_obj, info_in)
        self.snpeff_cat = search_result.group(1) 

    def return_attributes(self):
        return('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.chrom, self.pos, self.ref, self.alt, self.gene, self.amino_acid, self.aa_pos, self.sift_score, self.sift_median, self.n_seq, self.sift_cat, self.snpeff_cat))


def parse_sift(infile_in, snp_dict_in):
    '''
    function that loops through sift-annotated 
    vcf and parses each line, creating new
    AnnotatedSNP objects and adding them
    to the provided dictionary.
    '''
    for line in infile_in:
        if line[0] != "#":
            line = line.split("\t")
            chrom = line[0]
            site = line[1]
            ref = line[3]
            alt = line[4]
            info = line[7]
            #import pdb; pdb.set_trace()
            new_snp = AnnotatedSNP(chrom, site, ref, alt, info)
            snp_dict_in[site] = new_snp
    return snp_dict_in

def parse_snpeff(infile_in, snp_dict_in):
    '''
    function that loops through snpeff-annotated 
    vcf and parses each line, creating new
    AnnotatedSNP objects and adding them
    to the provided dictionary.
    '''
    for line in infile_in:
        if line[0] != "#":
            line = line.split("\t")
            site = line[1]
            info = line[7]
            snp_dict_in[site].parse_snpeff_info(info)
    return snp_dict_in

def write_outfile(snp_dict_in, outfile_in):
    outfile.write('CHROM\tPOS\tREF\tALT\tGENE\tAA\tAAPOS\tSIFTSCORE\tSIFTMED\tNSEQ\tSIFTCAT\tSNPEFFCAT\n')
    for snp in snp_dict_in:
        var_line = snp_dict_in[snp].return_attributes()
        outfile_in.write(var_line)
        

if __name__ == '__main__':
    # parse sift vcf
    print("sys.argv length: " + str(len(sys.argv)))
    snp_dict = {} 
    with open(sys.argv[1], 'r') as sift_infile:
        snp_dict = parse_sift(sift_infile, snp_dict)
    print("finished step 1\n")

    # parse snpeff vcf (if there is one)
    if len(sys.argv) > 3:
        with open(sys.argv[2], 'r') as snpeff_infile:
            snp_dict = parse_snpeff(snpeff_infile, snp_dict)
        print("\n\nfinished step 2\n")

    # print outfile
    with open(sys.argv[-1], 'w') as outfile:
        write_outfile(snp_dict, outfile)
