#!/usr/bin/env python

'''
Author: Stephen A. Sefick
Date: 20170314
Language: Python 2.7
Usage:
vcf2table_to_plot.py SNP.vcf SNP_table.out.txt
This script will take a SNP only vcf and strip out a number of annotations 
that may be of interest to plot and decide on hard filtering
'''
###########################################################################
##imports
import sys
import csv
import vcf
import re
###########################################################################

###########################################################################
##commands
##open both input and output files
##this automatically closes the files when the script is done
##with open('test3.vcf', 'r') as infile, open('test_out_py', 'w') as out_csv:
infile=sys.argv[1]
outfile=sys.argv[2]
with open(infile, 'r') as infile, open(outfile, 'w') as out_csv:    
    vcf_reader = vcf.Reader(infile)
    out = csv.writer(out_csv, delimiter='\t')
    header=['sample', 'CHROM', 'POS', 'ReadPosRankSum', 'MQRankSum', 'SOR', 'MQ', 'QD', 'GQ', 'DP', 'FS', 'DP_multi_sample', 'RGQ']
    out.writerow(header)
    for i in vcf_reader:
        for j in i.samples:
            if 'ReadPosRankSum' in i.INFO:
                RPRS=str(i.INFO['ReadPosRankSum'])
            else:
                RPRS='None'
            if 'MQRankSum' in i.INFO:
                MQRS=str(i.INFO['MQRankSum'])
            else:
                MQRS='None'
            ## if 'GQ' in str(j):
            ##     GQ=str(j['GQ'])
            if re.search(r'\bGQ', str(j)) is not None:
                GQ=str(j['GQ'])
            else:
                GQ='None'
            if 'DP' in str(j):
                DP=str(j['DP'])
            else:
                DP='None'
            if re.search(r'\bRGQ', str(j)) is not None:
                RGQ=str(j['RGQ'])
            else:
                RGQ='None'
            if 'SOR' in i.INFO:
                SOR=str(i.INFO['SOR'])
            else:
                SOR='None'
            if 'MQ' in i.INFO:
                MQ=str(i.INFO['MQ'])
            else:
                MQ='None'
            if 'QD' in i.INFO:
                QD=str(i.INFO['QD'])
            else:
                QD='None'
            if 'FS' in i.INFO:
                FS=str(i.INFO['FS'])
            else:
                FS='None'
            if 'DP' in i.INFO:
                DP_multi_sample=str(i.INFO['DP'])
            else:
                DP_multi_sample='None'    
            GQ_DP_SAMPLE=[j.sample, i.CHROM, i.POS, RPRS, MQRS, SOR, MQ, QD, GQ, DP, FS, DP_multi_sample, RGQ]
            out.writerow(GQ_DP_SAMPLE)            

###########################################################################
