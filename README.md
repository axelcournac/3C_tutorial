# 3C Analysis Tools

![quizzler_workflow](https://github.com/axelcournac/3C_analysis_tools.git/3C_analysis_tools/Capture du 2016-05-12 15:30:50.png)


This repository contains several codes for the processing, vizualisation and analysis of 3C/Hi-C data.
These codes were used during the INSERM workshop "Capturing chromosone conformation: toward a 3D view of genome regulation", May 9-13, Paris at Pasteur Institut.

For queries or help getting these running, you can contact me on mail or open an issue at the github repository.

## Dependencies

Scripts will run on OS X and other Unix-based systems. External dependencies should be installed somewhere on your `$PATH`.
It basically requires to have Python installed on your machine which commonly installed on Unix-based systems. 
For windows, you can have a look to https://www.python.org/downloads/windows/.

### Python 
* Numpy
* Matplotlib
* Scipy
* Biopython



### External programs

* `SRA tool` / [SRA](http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software)
* `Bowtie2 ` / [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)



## Session 1: Raw data extraction and alignment
Raw data are deposited on Short Read Archive server at the following address **http://www.ncbi.nlm.nih.gov/sra**.
We will take as example in this tutorial the 


We used an SRA executable called fastq-dump from SRA to extract and split both mates of a library.
```bash
/fastq-dump SRR639047 --split-3 -O /path_to_a_directory
```

# Alignment
We aligned the reads with iterative alignment procedure (scripts used: iterative_mapping.py and mapping.py of [hiclib] (https://bitbucket.org/mirnylab/hiclib) ). We modified the script mapping.py to add a value threshold on the mapping quality (30 or 40). The modification is given in the python directory. 


Before aligning the reads on a referecne genome, you need to index it:

```bash
bowtie2-build chr1.fa, chr2.fa name_index
```
You can also download the most used ones on the following link: **http://bowtie-bio.sourceforge.net/bowtie2/index.shtml**

Examples of lines used to launch the alignment procedure:
```bash
bowtie2 -x indices_genomes/sacCer3/sacCer3 -p6 --sam-no-hd --sam-no-sq --quiet --local --very-sensitive-local \
-S p1.sam /media/03b8b079-9d7a-4162-8201-6dd5d9923f62/2013/11_05_2013_Hi_Seq_MM/sequencage_nov2013/RSG6_L6/seq/BC76_CTGT.dat.end1.pcrfree

bowtie2 -x indices_genomes/sacCer3/sacCer3 -p6 --sam-no-hd --sam-no-sq --quiet --local --very-sensitive-local \
-S p2.sam /media/03b8b079-9d7a-4162-8201-6dd5d9923f62/2013/11_05_2013_Hi_Seq_MM/sequencage_nov2013/RSG6_L6/seq/BC76_CTGT.dat.end2.pcrfree
```


We convert the output of the aligment which is in HDF5 format into text file with the script convert_HDF5_txt.bh:

```bash
bash convert_HDF5_txt.bh /path_of_the_bank_of_the_output_of_aligment
```

## Filtering of the data


We then removed PCR duplicates. This is done using the C code pcr_duplicate.c (using hash table for C with the library uthash-1.9.6).
```bash
gcc pcr_duplicate.c
./a.out
```



## Normalization of the data
We used the normalization procedure called SCN (presented in http://www.biomedcentral.com/1471-2164/13/436) that we implemented in C. 
With dynamic allocation of memory, C language allows us to allocate big matrices (100000 x 100000) in a station with 50G of ram. The normalization is done in the code colocalization_cover.c 


This program takes several arguments in input: 
- a file of output of the aligment 
- a file containing the information concerning a biological parameter of the bins used for the analysis, i.e could be the GC content, the reads coverage, the chromosomes distribution that will be used to make random sets with a particular null model. 
- the thresholds [min - max] used for the normalization procedure that will exclude all bins whose norm is outside the ranges. It will notably filters poor interacting bins. 

To calculate the bins coverage from a file of output of alignment, we use the C code bins_coverage.c.
This code takes as input a file that contains coordinates of the bins and an alignment output file. 

To have bins that are enriched with a certain repeat (repeat_masker_2013.dat22 file), we used the C code binnage_distrib.c which calculate the p-value of the null model of distribution (binomial law) associated with every bins. The library libRmath must be installed to use R functions in C. 
To compile : 
```bash
gcc binnage_distrib.c -lRmath  -lm
```



### Scripts used for Figures 3





