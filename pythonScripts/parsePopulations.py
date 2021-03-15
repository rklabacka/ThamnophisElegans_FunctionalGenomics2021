#! usr/bin/env/ python
''' Script to parse individuals into
files for respective populations '''

#### Imports
import sys
import re
############

pop_set = set()
re_string = r'(?:MAR)|(?:ELF)|(?:PIK)|(?:MAH)|(?:PAP)|(?:NAM)|(?:MER)|(?:STO)|(?:SUM)|(?:PVM)|(?:RON)|(?:CHR)|(?:ROC)'
re_obj = re.compile(re_string)

def pairwise_populations(populationSet):
    for counter_a, pop_a in enumerate(populationSet):
#        print("pop_a: " + pop_a + ", counter_a: " + str(counter_a))
        for counter_b, pop_b in enumerate(populationSet):
#            print("pop_b: " + pop_b + ", counter_b: " + str(counter_b))
            if counter_b < counter_a:
                with open(pop_a + ".txt", 'r') as pop_a_file:
                    with open(pop_b + ".txt", 'r') as pop_b_file:
                        with open(pop_a + "_" + pop_b + ".txt", 'w') as pairwise_pops:
                            pairwise_pops.write(pop_a_file.read())
                            pairwise_pops.write(pop_b_file.read())

                
def all_populations():
    with open(sys.argv[1], 'r') as in_stream:
        for line in in_stream:
            re_search = re_obj.search(line)
            pop_set.add(re_search.group())
            with open(re_search.group() + ".txt", 'a') as out_stream:
                out_stream.write(line)

def combine_two_pops(pop1, pop2):
    with open(pop1 + ".txt", 'r') as pop1_file:
        with open(pop2 + ".txt", 'a') as pop2_file:
            pop2_file.write(pop1_file.read())

if __name__ == '__main__':
    all_populations()
    combine_two_pops("PIK", "MAR")
    pairwise_populations(pop_set)
