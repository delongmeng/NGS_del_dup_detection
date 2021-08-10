
# Genetic deletion and duplication calling

## Background

A disorder called Deldupemia is caused by the CNSL gene. High sequence homology in the CNSL gene region results in common deletions and duplications. The deletion/duplication breakpoints can vary from sample to sample, and it has been hypothesized that the breakpoints correspond with ethnicity. This “cnsl_data.csv.gz” dataset contains NGS reads at 100 different hybrid capture probes of 10,000 samples spanning different ethnicities. Of the probes, 50 probes are in the CNSL region, and 50 outside of the region. It's known that there are 4 breakpoint positions, and they correspond to the 32-38, 27-34, 20-40, and 10-40 probes in the CNSL gene region. 

This program processes the original reads, and annotate the genotype of CNSL gene for each sample, and save the resulted file. In addtion, it will also summarize the annotations by ethnicity for further analysis. 

## Run the program

- Environment  
Python, Pandas, Numpy
- Run the code (by default, this will take `cnsl_data.csv.gz` as input file):  
`python del_dup_analysis.py`
- Alternatively, you can specify the path of an input file of similar structure:    
`python del_dup_analysis.py --input_path xxxxxx`

## Approach

Please refer to the jupyter notebook file `CNSL_del_dup_analysis.ipynb` for step by step explanation of the process and visualization of the data. 

Briefly, I performed sample normalization first, based on the average reads of the 50 probes within the nonCNSL region, to adjust for the variability of extraction efficiency in the lab and error in the quantification of DNA libraries. For each sample, the reads will be divided by the average of the reads for all of the probes on the nonCNSL gene region. The assumption here is the average reads of the 50 nonCNSL region (they serve as control) should be the same from sample to sample, if there were no error in the extraction process.   

Next, I performed probe normalization because there is huge variation from probe to probe for a certain sample (due to different efficiency in DNA capture). The assumption here is, if there were no difference in the DNA capture efficiency of each probe, the average reads across all samples for each probe should be similar. Note that the deletion and duplication will for sure affect the reads of the probes within the region. We can ignore this effect because: 1) the overall deletion or duplication frequency is not high; 2) deletion and duplication will have opposite effects on the signal and thus cancel out each other.  

After normalization, 1 indicates CN = 2, 0.5 indicates a deletion (CN = 1), and 1.5 indicates a duplication (CN = 3). By data visualization, I notice that some of the probes have very high variation, which is not expected. These problematic probes don't form any contiguous stretch of at least four probes, so I selectively ignored them during the variation calling process.

For variation annotation, I check for each sample the normalized reads of the 50 CNSL gene probes. For each potential breakpoint (check from the biggest to the smallest range) I calculate the average reads and determine if it's a duplication (around 1.5) or deletion (around 0.5). See below an example of calling a deletion of breakpoint probe_27-34 (sample #163).

<img src ='/var/folders/tb/gs8f8b8n7mx6zhvfb0pftvqm0000gn/T/com.evernote.Evernote/WebKitDnD.h6MHNn/DC6F60D8-43EE-4BF9-9C9B-CA974F18A146.png' alt='image1' width='800'>



