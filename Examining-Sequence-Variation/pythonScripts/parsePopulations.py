#! usr/bin/env/ python
''' Script to parse individuals into
files for respective populations '''

#### Imports
import sys
import re
############


def all_populations(in_stream, re_obj):
    "create text file for each population"
    for line in in_stream:
        re_search = re_obj.search(line)
        pop_set.add(re_search.group())
        with open(re_search.group() + ".txt", 'a') as out_stream:
            out_stream.write(line)

def pairwise_populations(populationSet):
    for counter_a, pop_a in enumerate(populationSet):
#        print("pop_a: " + pop_a + ", counter_a: " + str(counter_a))
        for counter_b, pop_b in enumerate(populationSet):
#            print("pop_b: " + pop_b + ", counter_b: " + str(counter_b))
            if counter_b < counter_a:
                with open(pop_a + ".txt", 'r') as pop_a_file:
                    with open(pop_b + ".txt", 'r') as pop_b_file:
                        with open(pop_a + "_" + pop_b + ".txt", 'w') as pairwise_pops:
                            with open(pop_a + "_" + pop_b + ".pix", 'a') as pairwise_pix:
                                for line in pop_a_file:
                                    pairwise_pops.write(line)
                                    pairwise_pix.write(line.strip() + "\t" + pop_a + "\n")
                                for line in pop_b_file:
                                    pairwise_pops.write(line)
                                    pairwise_pix.write(line.strip() + "\t" + pop_b + "\n")
                
def combine_two_pops(pop1, pop2):
    with open(pop1 + ".txt", 'r') as pop1_file:
        with open(pop2 + ".txt", 'a') as pop2_file:
            pop2_file.write(pop1_file.read())

if __name__ == '__main__':
    pop_set = set()
    # Establish regular expression pattern for populations
    re_string = r'(?:MAR)|(?:ELF)|(?:PIK)|(?:MAH)|(?:PAP)|(?:NAM)|(?:MER)|(?:STO)|(?:SUM)|(?:PVM)|(?:RON)|(?:CHR)|(?:ROC)'
    # Compile regular expression pattern into regex object
    re_obj = re.compile(re_string)
    # Create text file for each population
    with open(sys.argv[1], 'r') as in_stream:
        all_populations(in_stream, re_obj)
    # Combine pik and mar files
    combine_two_pops("PIK", "MAR")
    # Create text file for each pairwise population comparison
    pairwise_populations(pop_set)
