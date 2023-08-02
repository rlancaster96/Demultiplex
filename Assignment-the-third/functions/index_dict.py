#!/usr/bin/env python

import bioinfo

index_nam = {"B1":"GTAGCGTA",
"A5":"CGATCGAT",
"C1":"GATCAAGG",
"B9":"AACAGCGA",
"C9":"TAGCCATG",
"C3":"CGGTAATC",
"B3":"CTCTGGAT",
"C4":"TACCGGAT",
"A11":"CTAGCTCA",
"C7":"CACTTCAC",
"B2":"GCTACTCT",
"A1":"ACGATCAG",
"B7":"TATGGCAC",
"A3":"TGTTCCGT",
"B4":"GTCCTAAG",
"A12":"TCGACAAG",
"C10":"TCTTCGAC",
"A2":"ATCATGCG",
"C2":"ATCGTGGT",
"A10":"TCGAGAGT",
"B8":"TCGGATTC",
"A7":"GATCTTGC",
"B10":"AGAGTCCA",
"A8":"AGGATAGC"}

def makedict(index_names: dict) -> dict:
    '''Takes a dictionary with index strings as values and returns a new dictionary with keys = index string, values = reverse compliment of index string.'''
    index_dict={}
    for index in index_names:
        indexrev = bioinfo.rev_comp(index_names.get(index))
        index_dict[index_names.get(index)] = indexrev
    return index_dict

print(bioinfo.makedict(index_nam))