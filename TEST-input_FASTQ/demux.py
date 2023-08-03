#!/usr/bin/env python 

#import modules
import itertools
import bioinfo
import gzip
import argparse
import itertools


#set up argparse 
def get_args():
    parser = argparse.ArgumentParser(description="Demultiplexes samples")
    # parser.add_argument("-i", "--indexfile", help="Tab separated index file", required=True, type=str)
    parser.add_argument("-R1", "--read1", help="Specify the file name for read 1", required=True, type=str)
    parser.add_argument("-R2", "--index1", help="Specify the file name for index 1", required=True, type=str)
    parser.add_argument("-R3", "--index2", help="Specify the file name for index 2", required=True, type=str)
    parser.add_argument("-R4", "--read2", help="Specify the file name for read 2", required=True, type=str)
    parser.add_argument("-i", "--indfile", help="Specify the file containing indexes", required=True, type=str)
    parser.add_argument("-c", "--cutoff", help="Specify quality cutoff for indexes", required=False, type=int)
    return parser.parse_args()

args = get_args()

#read indexes file, make dictionary of that data 
indexfiledict: dict = {}

with open(args.indfile, "r") as fh:
    fh.readline() #read the header line first to not include it
    for line in fh: #then for every other line
        line = line.strip()
        line = line.split()
        indexfiledict[line[4]] = line #make the indexfiledict: key is index seq, value is all data on line in a list

#Make dictionaries: one containing names and indexes, and one containing reverse compliment of indexes
index_names: dict = {}
for a in indexfiledict:
    index_names[indexfiledict.get(a)[0]] = indexfiledict.get(a)[4] #make a dictionary with sample : index

# index_names: dict = {"A1":"ATCG", "B2":"GGGG"} #unit test
rcindex: dict = bioinfo.makedict(index_names) #Takes a dictionary with index strings as values and returns a new dictionary with keys = index string, values = reverse compliment of index string.

#make dictionary for counting how many times we encounter an index (for percent of samples)
indexcounters: dict = {}
for a in index_names:
    indexcounters[index_names[a]]=0
totalsamples = 0 #initialize counter for total samples encountered

#make dictionary for counting how many times we encounter a PAIR of indexes (overall index swapping)
iterlist = []
for a in index_names:
    iterlist.append(index_names[a])
combos = list(itertools.product(iterlist, repeat=2))

paircounters: dict = {}
for a in combos:
    paircounters[a] = 0



# open read files 
read1 = open(args.read1, "r")
read2 = open(args.read2, "r")
index1 = open(args.index1, "r")
index2 = open(args.index2, "r")

# open the 4 static write files
R1_hopped = open("R1_hopped.fastq", "w")
R4_hopped = open("R4_hopped.fastq", "w")
R1_unknown = open("R1_unknown.fastq", "w")
R4_unknown = open("R4_unknown.fastq", "w")

handles: list = [R1_hopped, R4_hopped, R1_unknown, R4_unknown] #make a list of all filehandles to close at the end 
indexhandledict: dict = {}
#make filehandles and dictionary of filehandles for paired indexes where key = index, value = filehandle
for index in index_names:
    filehandle1 = "fh_R1_" + index_names[index]
    filename1 = "R1_index_" + index_names[index] + ".fastq"
    
    filehandle2 = "fh_R4_" + index_names[index]
    filename2 = "R4_index_" + index_names[index] + ".fastq"
    
    filehandle1 = open(filename1, "w")
    filehandle2 = open(filename2, "w")

    handles.append(filehandle1) #append filehandle names to the handles file for final closing
    handles.append(filehandle2)

    pack: tuple = (filehandle1, filehandle2) #add key = index, filehandles1 and 2 = values for each index. use this for writing paired indexes later.
    indexhandledict[index_names.get(index)] = pack

#initialize counters for matching, hopped, and unknown indexes
matching = 0
hopped = 0
unknown = 0

#start looping through files
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
    rcind2 = bioinfo.rev_comp(ind2) #make reverse complement of index 2
    appender = " " + ind1 + "-" + rcind2
    read1_header = read1_header + appender
    read2_header = read2_header + appender
    
    index1_qscore = index1_qscore.strip()
    index2_qscore = index2_qscore.strip()

    #1 : QSCORE FILTERING 
    #*************************************
    #get lowest qscore from indexes 
    lowscore1 = 100
    lowscore2 = 100
    
    for a in index1_qscore:
        qscore1 = bioinfo.convert_phred(a)
        if qscore1 < lowscore1:
            lowscore1 = qscore1
    for a in index2_qscore:
        qscore2 = bioinfo.convert_phred(a)
        if qscore2 < lowscore2:
            lowscore2 = qscore2
    # write to unknown file if it doesn't make the cutoff
    if lowscore1 < args.cutoff or lowscore2 < args.cutoff: 
        unknown +=1
        R1_unknown.write(f"{read1_header}\n{read1_seq}{read1_comment}{read1_qscore}")
        R4_unknown.write(f"{read2_header}\n{read2_seq}{read2_comment}{read2_qscore}")
    
    #2 : DICTIONARY COMPARISON
    #*************************************
    #write to unknown file if index1 or rev comp of index 2 is not a key in the rcindex dictionary:
    elif ind1 not in rcindex or rcind2 not in rcindex:
        unknown +=1
        R1_unknown.write(f"{read1_header}\n{read1_seq}{read1_comment}{read1_qscore}")
        R4_unknown.write(f"{read2_header}\n{read2_seq}{read2_comment}{read2_qscore}")
    
    #3 : PAIRED INDEXES
    #*************************************
    #write to paired indexes file if indexes match 
    elif ind1 == rcind2:
        temptuple = (ind1, rcind2)
        paircounters[temptuple] +=1 #add count of 1 to dictionary counting pairs
        matching +=1 #add count of 1 to count of matched 
        indexcounters[ind1] += 1 #add count of 1 to individual index counter
        totalsamples +=1
        #get filehandles tuple
        fh1, fh2 = indexhandledict[ind1]
        fh1.write(f"{read1_header}\n{read1_seq}{read1_comment}{read1_qscore}")
        fh2.write(f"{read2_header}\n{read2_seq}{read2_comment}{read2_qscore}")

    #4 : HOPPED FILES
    #*************************************
    else:
        hopped +=1
        temptuple = (ind1, rcind2)
        try:
            paircounters[temptuple] +=1 
        except:
            print(hopped, "NOT FOUND!",temptuple)
        R1_hopped.write(f"{read1_header}\n{read1_seq}{read1_comment}{read1_qscore}")
        R4_hopped.write(f"{read2_header}\n{read2_seq}{read2_comment}{read2_qscore}")
file=open("LESLIE","w")
for thing in paircounters:
    print(thing, file=file)

#close files
for handle in handles:
    handle.close()

# Markdown output file
#*********************

#Calculate percent of reads from each sample:
percentsample: dict = {} #make a dictionary with index=key, value = tuple(percent, sample#)
for a in indexcounters:
    percent = ((indexcounters[a])/totalsamples)*100
    percentsample[a] = [percent, indexfiledict[a][0]]

with open("Demux_Summary.md", "w") as wh:
    wh.write(f"Demultiplexing Summary\n")
    wh.write(f"Hopped: {hopped}\nMatching: {matching}\nUnknown: {unknown}")
    wh.write(f"\n\nPercent of reads from each sample\n")
    wh.write(f"Sample\tRead index\tPercent\n")
    for a in percentsample:
        wh.write(f"{percentsample[a][1]}\t{a}\t{percentsample[a][0]}\n")
    wh.write(f"\nOverall amount of index swapping\n")
    wh.write(f"Pair combination\tFrequency\n")
    for a in paircounters:
        wh.write(f"{a}\t{paircounters[a]}\n")


#testing delete later 
# print(index_names)
# print(indexcounters)
# print(len(index_names))
# print(len(indexcounters))

# print(paircounters)
# print(indexhandledict)