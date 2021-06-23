#! usr/bin/env/ python
import sys

Genes_in = open(sys.argv[1], 'r')
GFF_in = open(sys.argv[2], 'r')
Genes_var = []
GFF_out = open(sys.argv[3], 'w')
VariableGenes_out = open(sys.argv[4], 'w')

for line in Genes_in:
  line = line.split("\t")
  if int(line[1]) > 0:
    Genes_var.append(line[0])
    VariableGenes_out.write(line[0] + "\n")

for line in GFF_in:
  gene = line.split("gene=")[1]
  gene = gene.split(";")[0]
  if gene in Genes_var:
    GFF_out.write(line)


Genes_in.close()
GFF_in.close()
GFF_out.close()
VariableGenes_out.close()
