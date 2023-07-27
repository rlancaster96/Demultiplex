#!/usr/bin/env python

#import modules
import bioinfo
import matplotlib.pyplot as plt
import argparse
import gzip

#set up args for input from command line
def get_args():
    parser = argparse.ArgumentParser(description="Takes a fastq file, read length, and output filename and creates a histogram of mean distribution of quality scores at each base pair position")
    parser.add_argument("-l", "--readlen", help="Specify the read length", required=True, type=int)
    parser.add_argument("-f", "--filename", help="Specify the file name", required=True, type=str)
    parser.add_argument("-o", "--writename", help="Specify output file name", required=False, type=str)
    
    return parser.parse_args()

args = get_args()
args.readlen = args.readlen + 1 #read length plus one for exclusive range

#Modify init_list from PS4 to make length of the list variable (accounts for indexes too). 
def init_list(lst: list, rangemax: int, value: float=0.0) -> list:
    '''This function takes an empty list, an exlusive top end of the range value,
    and will populate it with the value passed in "value". If no value is passed, 
    initializes list with values of 0.0.'''
    for i in range(1,rangemax):
        lst.append(value)
    return lst

my_list=[]
init_list(my_list, args.readlen)

#Modify populate_list from PS4
def populate_list(file: str, rangemax: int) -> tuple[list, int]:
    '''Opens the gzipped file passed, converts phred quality scores from letters to numbers, and returns a list of sums of all quality scores for each base pair position.'''
     # initialize the empty list
    with gzip.open(file, "rt") as fh: #read file, uncomment this for gzipped files
    # with open(file, "rt") as fh: #read file, uncomment this for regular files
        i = 0 # initialize counter
        for line in fh:
            
            i+=1 #step counter for each line
            line = line.strip('\n') #strip newline char
            if i%4 == 0: #every fourth line (the quality score line)
                for a in range(len(line)): #for every character on the line
                    my_list[a] += bioinfo.convert_phred(line[a]) #the value of quality score at each base pair position is added
    return(my_list,i)

my_list, num_lines = populate_list(args.filename, args.readlen)

#make list of means
for i in range (len(my_list)): #for each value in the list
    my_list[i] = (my_list[i]/(num_lines/4)) # replace that value with the value divided by number of quality score lines in the file (num_lines/4)

#optional print step: uncomment below to print distribution to terminal.
print(f"#Base pair\tMean Quality Score")
for i in range (len(my_list)):
    print(f"{i+1}\t{my_list[i]}") #print in the requested format

#plot hist
plt.plot(my_list)
plt.title("Average Quality Score per Base Pair Position in the Read")
plt.xlabel("Base Pair #")
plt.ylabel("Average Quality Score")
plt.savefig(f"{args.writename}.png")