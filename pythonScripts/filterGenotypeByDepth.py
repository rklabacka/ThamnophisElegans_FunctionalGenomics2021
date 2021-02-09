#! usr/bin/env python
import sys
import vcf
import os

vcf_reader = vcf.Reader(open(sys.argv[1], 'r'))
vcf_writer = vcf.Writer(open(sys.argv[2], 'w'), vcf_reader)
for record in vcf_reader:
    if record.INFO['DP']>20:
        vcf_writer.write_record(record)
