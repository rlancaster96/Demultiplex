# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

https://github.com/rlancaster96/Demultiplex/blob/master/Assignment-the-first/qscore_dist.py

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read 1 | 101 | Phred +33 |
| 1294_S1_L008_R2_001.fastq.gz | index 1 | 8 | Phred +33 |
| 1294_S1_L008_R3_001.fastq.gz | index 2 | 8 | Phred +33 |
| 1294_S1_L008_R4_001.fastq.gz | read 2 | 101 | Phred +33 |

2. Per-base NT distribution

   i. Histograms:


   
![R1_hist](https://github.com/rlancaster96/Demultiplex/assets/136844363/82f72a0a-cfb2-4edf-83ee-23431fa43d51)

![R2_hist](https://github.com/rlancaster96/Demultiplex/assets/136844363/df00c274-20ef-4ea5-9c7c-b78ce9d89b84)

![R3_hist](https://github.com/rlancaster96/Demultiplex/assets/136844363/6b12e9ae-ecd2-41cf-ad2d-684199663e7b)

![R4_hist](https://github.com/rlancaster96/Demultiplex/assets/136844363/6b4cf3ee-3a66-4e9a-8b88-f36f6e745e94)



ii. I am basing my cutoffs on the means we just calculated. I would use a quality score cutoff of 33 (binned 30-34) for *any single base pair position* for my indexes and a cutoff of 27 (binned 25-29) for the *mean quality score across the read* of my biological reads (since they haven't yet been trimmed, I'd go a little lower than 33). The lowest mean quality score per base was between 30 and 32 for both the indexes and my biological reads. However, I want to have a stricter quality filter for indexes since I want to avoid index-swapping issues. These issues have worse ramifications for downstream analysis (including unrelated genomic data in my analysis) than for a few bad-quality nucleotides in a bioread. 

iii. 3976613 in R2, and 3328051 for R3. A total of 7304664. 
	    ```
	    zcat <filename> | sed -n '2~4p' | grep -E "[N]+" | wc -l 
	    ```

## Part 2
1. Define the problem

This library was prepared with dual paired indexes. I need to de-multiplex our sequence files and determine the amount of index hopping. We have to filter sequencing data for quality, separating out passing paired indexes from indexes that are too low quality, don't match the expected indexes, or have hopped/swapped.
Sequence runs come in 4 files: the first read (R1), the first index (R2), the second index (R3), and the second read (R4). I have to parse these four files at the same time to compare the first and second indices - if one index is the reverse compliment of the other, matches an index in the provided index list, and is of sufficient quality (doesn't contain "N"s), then it is a passing index pair and read1/read2 can be written to my passing files for that specific index pair. Indexes that are hopped go into the r1/r2 hopped files, and those not present in the provided index list or are of insufficient quality also go into the read1/read2 low quality files.

2. Describe output

The outputs are FASTQ files. The headers for the records should have the index pair appended to them (in the format of "index1"-"reverse compliment of index2" for ease of reading). There is one FASTQ file per index per read for correctly paired indexes that pass the quality filter. Swapped indexes will be appended to either the read 1 swapped or read 2 swapped file, and low quality indexes will be appeneded to either the read 1 low-q or read 2 low-q files. The total number of outputs given example indexes *n* and *x* will be:
- read 1 index *n*
- read 2 index *n*
- read 1 index *x*
- read 2 index *x*
- read 1 - all swapped
- read 2 - all swapped
- read 1 - all low-quality
- read 2 - all low-quality 

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).

4. Pseudocode: 
```
make dict of indexes and their reverse compliments 
indexdict = {index1 : RC-index1, index2 : RC-index2, ...}

initialize counters for matching, hopped, and unknown indexes

open R1, R2, R3, R4 as filehandles
store header from R1 and R4 in temporary memory 
store rest of record from R1 and R4 in temporary memory 
store index from R2 and R3 in temporary memory as variables
index 1 = R2 index
index 2 = reverse compliment of R3 index 
append index1-index2 to headers 
use convert_phred function to convert string to quality score and find lowest quality score of index 1 and index 2.

begin quality filtering: 
if index1 doesnt match to dictionary key (known index) OR index2 doesn't match to dictionary value (reverse compliment of known index): 
	write modified header and rest of record to the R1 and R2 lowquality files
	unknown indexes += 1
elif index1 OR index2 lowest quality score is < cutoff: 
	write modified header and rest of record to the R1 and R2 lowquality files
	unknown indexes += 1
elif index1 == index2:
	write modified header and rest of record to the R1 and R2 passing files for the specific index
	matching indexes += 1
else: 
	write modified header and rest of record to the R1 and R2 hopped files 
	hopped indexes += 1

return values for unknown, matching, and hopped indexes
```

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement

I have four high-level functions I'm planning on using. Two are new and two are already in my bioinfo.py module.

```
def rev_comp(nucleotides: string) -> string: 
'''Takes a nucleotide string and returns the reverse compliment of that nucleotide string'''
return reversecompliment
Input: ACTG
Output: CAGT
```

```
def index_dict(indexlist: list) -> dict: 
'''Takes a list of indexes and returns a dictionary of key = "index" and value = "reverse compliment of index"'''
return index_dict
Input: [ACTG, AAAA, TTTT, CCCC]
Output: {ACTG : CAGT, AAAA : TTTT, TTTT : AAAA, CCCC : GGGG}
```

Already in bioinfo.py module: 

```
def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score. Uses a -33 phred convert.'''
    phred = int(ord(letter))-33  #Convert the letter to Unicode, convert to an integer, then subtract -33 to get the phred score
    return phred
Input: E
Output: 36
```

