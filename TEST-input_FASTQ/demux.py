#!/usr/bin/env python 

#import modules
import itertools
import bioinfo
import gzip
import argparse

#set up argparse 
def get_args():
    parser = argparse.ArgumentParser(description="Demultiplexes samples")
    # parser.add_argument("-i", "--indexfile", help="Tab separated index file", required=True, type=str)
    parser.add_argument("-R1", "--read1", help="Specify the file name for read 1", required=True, type=str)
    parser.add_argument("-R2", "--index1", help="Specify the file name for index 1", required=True, type=str)
    parser.add_argument("-R3", "--index2", help="Specify the file name for index 2", required=True, type=str)
    parser.add_argument("-R4", "--read2", help="Specify the file name for read 2", required=True, type=str)
    # parser.add_argument("-o", "--writename", help="Specify output file name", required=False, type=str)
    return parser.parse_args()

args = get_args()

#Make dictionaries: one containing names and indexes, and one containing reverse compliment of indexes
# index_names: dict = {"B1":"GTAGCGTA", "A5":"CGATCGAT", "C1":"GATCAAGG",
# "B9":"AACAGCGA", "C9":"TAGCCATG", "C3":"CGGTAATC", 
# "B3":"CTCTGGAT", "C4":"TACCGGAT", "A11":"CTAGCTCA",
# "C7":"CACTTCAC", "B2":"GCTACTCT", "A1":"ACGATCAG",
# "B7":"TATGGCAC", "A3":"TGTTCCGT", "B4":"GTCCTAAG",
# "A12":"TCGACAAG", "C10":"TCTTCGAC", "A2":"ATCATGCG",
# "C2":"ATCGTGGT", "A10":"TCGAGAGT", "B8":"TCGGATTC",
# "A7":"GATCTTGC", "B10":"AGAGTCCA", "A8":"AGGATAGC"}

index_names: dict = {"A1":"ATCG", "B2":"GGGG"}
rcindex: dict = bioinfo.makedict(index_names)

# open files 
read1 = open(args.read1, "r")
read2 = open(args.read2, "r")
index1 = open(args.index1, "r")
index2 = open(args.index2, "r")

while True: 

#store 4 records at a time (one from each file) to temporary memory 
    
    read1_header = read1.readline()
    if read1_header == "": #once we reach the end of the file (no more headers), break loop
        break
    read1_header = read1_header.strip()
    read1_seq = read1.readline()
    read1_comment = read1.readline()
    read1_qscore = read1.readline()

    read2_header = read2.readline()
    read2_header = read2_header.strip()
    read2_seq = read2.readline()
    read2_comment = read2.readline()
    read2_qscore = read2.readline()

    index1_header = index1.readline()
    index1_seq = index1.readline()
    index1_comment = index1.readline()
    index1_qscore = index1.readline()

    index2_header = index2.readline()
    index2_seq = index2.readline()
    index2_comment = index2.readline()
    index2_qscore = index2.readline()

    #strip newline characters and make index variables
    ind1 = index1_seq.strip()
    ind2 = index2_seq.strip()
    rcind2 = bioinfo.rev_comp(ind2)
    appender = " " + ind1 + "-" + rcind2
    read1_header = read1_header + appender
    read2_header = read2_header + appender
    
    index1_qscore = index1_qscore.strip()
    index2_qscore = index2_qscore.strip()


    #get lowest qscore from indexes 
    lowscore1 = 100
    lowscore2 = 100
    for a in index1_qscore:
        qscore1 = bioinfo.convert_phred(index1_qscore)
        if qscore1 < lowscore1:
            qscore1 = lowscore1
    for a in index2_qscore:
        qscore2 = bioinfo.convert_phred(index2_qscore)
        if qscore2 < lowscore2:
            qscore2 = lowscore2
    
    print(lowscore1, lowscore2)



#close files
read1.close()
read2.close()
index1.close()
index2.close()
