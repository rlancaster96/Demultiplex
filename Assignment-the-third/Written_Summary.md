**Demultiplexing Write-Up 8-11-2023**

For this assignment, I demultiplexed samples from a sequencing experiment performed by a previous BGMP cohort in 2017. 

I used a quality filter that used per-nucleotide quality filtering for the index sequences. The quality filter is adjustable using argparse, but for my final data I decided on a cutoff of 30. I checked that my program was returning expected values by running it with a cutoff of 2 exclusive (this is the same as only filtering out "N") and comparing my data to others who had used the same quality filtering parameters. 

Although I do think my cutoff of 30 per base is quite restrictive, I am retaining ~62% of reads for downstream analysis (~220 million total). The lowest quantity of reads I achieved for any sample was about 2 million (this was sample #7), and the highest about 50 million (sample #10). From preliminary searching, I found that this is sufficient for some downstream analysis, especially when using a panel like Illumina's TruSight RNA Pan Cancer, but not enough for something like de novo transcriptome assembly (recommended at least 100 million reads), and it does depend on the organism I'm looking at (source:https://web.genewiz.com/rna-seq-faq, https://knowledge.illumina.com/library-preparation/rna-library-prep/library-preparation-rna-library-prep-reference_material-list/000001243). If I proceeded to further analysis and found that my dataset was insufficient, I can easily change it to a lower quality cutoff to get more data. 

The weakness of my program is that it does not take advantage of the fact that we have 1 nucleotide of uncertainty that is allowable because of how our indexes are designed. I could retain even more data if I had compared indexes and allowed 1 nucleotide of variation between the observed sequence and the known index. I think this would be great to add as an option to pass the algorithm, especially if I am working with a low-quality, low-quantity dataset and want to retain as much as data as possible. 

**Observations of the Data**

Samples 10 (index #TACCGGAT) and 23 (index # TCTTCGAC) had the most representation in final matching reads. 

The index pairs with the most swapping events were these: 

```
('TATGGCAC', 'TGTTCCGT'), 58741 events
('TGTTCCGT', 'TATGGCAC'), 58569 events 
```


![Summary_Index_Swapping_Heatmap](https://github.com/rlancaster96/Demultiplex/assets/136844363/92f0472c-03bc-4513-93f6-c73b4cb9615f)


I had a total index hopping rate of ~0.15%, which is a typical rate and suggests that the amount of index hopping is normal. It is usually between 0.1-0.2% (source: https://www.illumina.com/techniques/sequencing/ngs-library-prep/multiplexing/index-hopping.html)
