#!/usr/bin/env python

import itertools

index_names: dict = {"B1":"GTAGCGTA", "A5":"CGATCGAT", "C1":"GATCAAGG",
"B9":"AACAGCGA", "C9":"TAGCCATG", "C3":"CGGTAATC", 
"B3":"CTCTGGAT", "C4":"TACCGGAT", "A11":"CTAGCTCA",
"C7":"CACTTCAC", "B2":"GCTACTCT", "A1":"ACGATCAG",
"B7":"TATGGCAC", "A3":"TGTTCCGT", "B4":"GTCCTAAG",
"A12":"TCGACAAG", "C10":"TCTTCGAC", "A2":"ATCATGCG",
"C2":"ATCGTGGT", "A10":"TCGAGAGT", "B8":"TCGGATTC",
"A7":"GATCTTGC", "B10":"AGAGTCCA", "A8":"AGGATAGC"}

iterlist = []
for a in index_names:
    iterlist.append(index_names[a])

combos = list(itertools.combinations_with_replacement(iterlist,2))
print(combos)