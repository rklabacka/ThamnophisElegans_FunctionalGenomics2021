#! usr/bin/env python
import sys

myList = []

with open(sys.argv[1], "r") as gff:
    haystack = gff.read()

if not haystack:
    sys.exit("Could not read haystack data :-(")

with open(sys.argv[2], "r") as genes:
    for needle in (line.strip() for line in genes):
        if needle in haystack:
            print(needle, ',found')
            if needle not in myList:
                myList.append(needle)
        else:
            print(needle, ',not found')

with open(sys.argv[3], "w") as outfile:
    for item in myList:
        outfile.write(item + "\n")
