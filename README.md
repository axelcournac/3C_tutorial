# 3C Analysis Tools

This repository contains scripts to recreate figures in paper **Mutiscales** (2016). For queries or help getting these running, you can contact me on mail or open an issue at the github repository.


## Dependencies

Scripts will run on OS X and other Unix-based systems. External dependencies should be installed somewhere on your `$PATH`.

### R packages

Lots of commonly-installed R packages are also used, including but not limited to: 

#### CRAN

* `calibrate`
* `verification`
* `ROCR`
* `libRmath`

#### Bioconductor

* `BSgenome` (and UCSC hg19)

### External programs

* `SRA tool` / [SRA](http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software)
* `Bowtie2 ` / [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* `ICE` / [hiclib](http://mirnylab.bitbucket.org/hiclib)
* `GNU Scientific Library` / [GSL](http://www.gnu.org/software/gsl/) 
* `uthash a hash table in C` / [uthash](https://troydhanson.github.io/uthash/userguide.html)


## Raw data
Raw data were deposited on Short Read Archive server at the following address, example for the first run of the bank 




We separate both ends of the reads using the command lines written in the script separate_mates.sh :
```bash
bash separate_mates.sh
```
We aligned the reads with iterative alignment procedure (scripts used: iterative_mapping.py and mapping.py of [hiclib] (https://bitbucket.org/mirnylab/hiclib) ). We modified the script mapping.py to add a value threshold on the mapping quality (30 or 40). The modification is given in the python directory. 
Example of lines used to launch the alignment procedure:
```bash
bank='/media/human/bank400260/';  echo $bank; 
fast1='/media/human/seq/SRR400260.fastq.end1'
fast2='/media/human/seq/SRR400260.fastq.end2' 
path_to_index='/home/axel/Bureau/python/fasta/human/hg19_python';
path_to_fasta='/home/axel/Bureau/python/fasta/human';
name_of_enzyme='HindIII'

python iterative_mapping.py  $bank $fast1 $fast2 $path_to_index $path_to_fasta $name_of_enzyme 
```

We convert the output of the aligment which is in HDF5 format into text file with the script convert_HDF5_txt.bh:

```bash
bash convert_HDF5_txt.bh /path_of_the_bank_of_the_output_of_aligment
```

## Filtering of the data



We then removed PCR duplicates i.e paires of reads that have exactly the same positions. This is done using the C code pcr_duplicate.c (using hash table for C with the library uthash-1.9.6).
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




## sessionInfo()

Below is an output of sessionInfo() for troubleshooting purposes, some loaded packages may not be required and likewise, some required packages may not be loaded. An exception caused by attached packages is likely due to version issues.
```r
R version 3.2.0 (2015-04-16)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

locale:
 [1] LC_CTYPE=fr_FR.utf8       LC_NUMERIC=C             
 [3] LC_TIME=fr_FR.utf8        LC_COLLATE=fr_FR.utf8    
 [5] LC_MONETARY=fr_FR.utf8    LC_MESSAGES=fr_FR.utf8   
 [7] LC_PAPER=fr_FR.utf8       LC_NAME=C                
 [9] LC_ADDRESS=C              LC_TELEPHONE=C           
[11] LC_MEASUREMENT=fr_FR.utf8 LC_IDENTIFICATION=C      

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base  
```

