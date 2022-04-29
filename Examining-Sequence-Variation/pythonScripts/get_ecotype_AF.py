#! usr/bin/env/ python
'''
Script for getting allele frequencies for each site of a VCF file
for ecotypes (each ecotype has an independent VCF file)

Input:
    Lakeshore.vcf (Ecotype 1 VCF)
    Meadow.vcf (Ecotype 2 VCF)

Output:
    EcotypicAlleleFrequencies.csv (comma-delimmited file with allele freqs)

Examples:
python get_ecotype_AF.py Lakeshore.vcf Meadow.vcf EcotypicAlleleFrequencies.csv
'''

import sys
import re
import argparse
import pysam

class VariableSite:
    def __init__(self, name_in, ref_in, alt_in):
        self.name = str(name_in)
        self.ref_allele = ref_in
        alt_in = re.sub(r'[()\',]','',alt_in)
        self.alt_allele = alt_in 
        self.lakeshore_af = 'NA'
        self.meadow_af = 'NA'

    def add_af(self, af, ecotype):
        af = str(af)
        af = re.sub(r'[(),]','',af)
        af = af[0:4]
        if ecotype == "Lakeshore":
            self.lakeshore_af = af
        elif ecotype == "Meadow":
            self.meadow_af = af
        else:
            print("ERROR: invalid ecotype")

def add_first_ecotype(site_dict, vcf_instream, ecotype):
    for line in vcf_instream:
        site_dict[line.pos] = VariableSite(line.pos, line.ref, str(line.alts))
        af = line.info['AF']
        add_allele_freq(site_dict, line.pos, af, ecotype)

def add_second_ecotype(site_dict, vcf_instream, ecotype):
    for line in vcf_instream:
        af = line.info['AF']
        add_allele_freq(site_dict, line.pos, af, ecotype) 

def add_allele_freq(site_dict, site, af, ecotype):
    site_dict[site].add_af(af, ecotype)

def write_outfile(site_dict, outstream):
    outstream.write("POS,REF,ALT,L_AF,M_AF")
    for site in site_dict:
        outstream.write(site_dict[site].name + ',' + 
            site_dict[site].ref_allele + ',' + 
            site_dict[site].alt_allele + ',' + 
            site_dict[site].lakeshore_af + ',' +
            site_dict[site].meadow_af + '\n')

def create_parser():
    parser = argparse.ArgumentParser(description = "Determines allele frequency of ecotype VCF files")
    parser.add_argument('-m', '--meadow', help = 'Required input VCF for meadow ecotype')
    parser.add_argument('-l', '--lakeshore', help = 'Required input VCF for lakeshore ecotype')
    parser.add_argument('-c', '--csv', help = 'Produce CSV with headers containing allele frequencies at each site for each ecotype')
    return parser.parse_args()

if __name__ == "__main__":
    assert len(sys.argv) == 7 
    args = create_parser()
    site_dict = {}
    with pysam.VariantFile(args.meadow, 'r') as meadow_in:
        add_first_ecotype(site_dict, meadow_in, "Meadow")
    with pysam.VariantFile(args.lakeshore, 'r') as lakeshore_in:
        add_second_ecotype(site_dict, lakeshore_in, "Lakeshore")
    with open(args.csv, 'w') as outfile:
        write_outfile(site_dict, outfile)
     
    

