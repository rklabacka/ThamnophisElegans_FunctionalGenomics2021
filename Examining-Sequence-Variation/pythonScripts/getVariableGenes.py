#! usr/bin/env/ python
import sys

Nvar_in = open(sys.argv[1], 'r')
Var_out = open(sys.argv[2], 'w')

leng = 0
past_gene = "start"

for line in Nvar_in:
    line = line.split("\t")
    if int(line[1]) > 0:
        print(str(line[1]))
        Var_out.write(line[0] + "\n")

Nvar_in.close()
Var_out.close()
