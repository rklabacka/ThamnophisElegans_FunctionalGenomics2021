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
        self.lakeshore_ac = 'NA'
        self.meadow_ac = 'NA'
        self.lakeshore_af = 'NA'
        self.meadow_af = 'NA'

    def add_af(self, ac, af, ecotype):
        "Add allele count and allele frequency values"
        af = self.modify_info_str(af)
        ac = self.modify_info_str(ac)
        if ecotype == "Lakeshore":
            self.lakeshore_ac = ac
            self.lakeshore_af = af
        elif ecotype == "Meadow":
            self.meadow_ac = ac
            self.meadow_af = af
        else:
            print("ERROR: invalid ecotype")

    def modify_info_str(self, x):
        "Modify AF or AC value from VCF using regex"
        x = str(x)
        x = re.sub(r'[(),]','',x)
        x = x[0:4]
        return x
        

def add_first_ecotype(site_dict, vcf_instream, ecotype):
    for line in vcf_instream:
        site_dict[line.pos] = VariableSite(line.pos, line.ref, str(line.alts))
        ac = line.info['AC']
        af = line.info['AF']
        add_allele_freq(site_dict, line.pos, ac, af, ecotype)

def add_second_ecotype(site_dict, vcf_instream, ecotype):
    for line in vcf_instream:
        ac = line.info['AC']
        af = line.info['AF']
        add_allele_freq(site_dict, line.pos, ac, af, ecotype) 

def add_allele_freq(site_dict, site, ac, af, ecotype):
    site_dict[site].add_af(ac, af, ecotype)

def write_outfiles(site_dict, outstream_all, outstream_abbrev):
    outstream_all.write("POS,REF,ALT,L_AC,M_AC,L_AF,M_AF\n")
    outstream_abbrev.write("POS,REF,ALT,L_AC,M_AC,L_AF,M_AF\n")
    for site in site_dict:
        if abs(float(site_dict[site].lakeshore_af) - float(site_dict[site].meadow_af)) > 0.3:
            write_lines_to_outfile(site_dict, site, outstream_abbrev)
        write_lines_to_outfile(site_dict, site, outstream_all)

def write_lines_to_outfile(site_dict, site, outstream):
    outstream.write(site_dict[site].name + ',' + 
        site_dict[site].ref_allele + ',' + 
        site_dict[site].alt_allele + ',' + 
        site_dict[site].lakeshore_ac + ',' +
        site_dict[site].meadow_ac + ',' +
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
    with open(args.csv, 'w') as outfile_all:
        with open("abbrev_" + args.csv, 'w') as outfile_abbrev:
            write_outfiles(site_dict, outfile_all, outfile_abbrev)
     
    

