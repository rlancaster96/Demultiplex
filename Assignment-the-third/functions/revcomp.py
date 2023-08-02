#!/usr/bin/env python

sequence = "TAAAAAAAANG"

def rev_comp(seq: str) -> str:
    '''Takes a dictionary called revdict (containing nucleotide complements) and a string. Returns reverse complement string. Revdict defined in bioinfo module.'''
    revlist = (list(reversed(seq)))
    revcomp = []
    for nuc in revlist:
        revcomp.append(revdict.get(nuc))
    revcomp = ''.join(revcomp)
    return revcomp
revdict = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
