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

### Alignment

Before aligning the reads on a reference genome, you need to index it:

```bash
bowtie2-build chr10.fa,chr11.fa,chr11_gl000202_random.fa,chr12.fa,chr13.fa,chr14.fa,chr15.fa,chr16.fa,chr17_ctg5_hap1.fa,chr17.fa,chr17_gl000203_random.fa,chr17_gl000204_random.fa,chr17_gl000205_random.fa,chr17_gl000206_random.fa,chr18.fa,chr18_gl000207_random.fa,chr19.fa,chr19_gl000208_random.fa,chr19_gl000209_random.fa,chr1.fa,chr1_gl000191_random.fa,chr1_gl000192_random.fa,chr20.fa,chr21.fa,chr21_gl000210_random.fa,chr22.fa,chr2.fa,chr3.fa,chr4_ctg9_hap1.fa,chr4.fa,chr4_gl000193_random.fa,chr4_gl000194_random.fa,chr5.fa,chr6_apd_hap1.fa,chr6_cox_hap2.fa,chr6_dbb_hap3.fa,chr6.fa,chr6_mann_hap4.fa,chr6_mcf_hap5.fa,chr6_qbl_hap6.fa,chr6_ssto_hap7.fa,chr7.fa,chr7_gl000195_random.fa,chr8.fa,chr8_gl000196_random.fa,chr8_gl000197_random.fa,chr9.fa,chr9_gl000198_random.fa,chr9_gl000199_random.fa,chr9_gl000200_random.fa,chr9_gl000201_random.fa,chrM.fa,chrUn_gl000211.fa,chrUn_gl000212.fa,chrUn_gl000213.fa,chrUn_gl000214.fa,chrUn_gl000215.fa,chrUn_gl000216.fa,chrUn_gl000217.fa,chrUn_gl000218.fa,chrUn_gl000219.fa,chrUn_gl000220.fa,chrUn_gl000221.fa,chrUn_gl000222.fa,chrUn_gl000223.fa,chrUn_gl000224.fa,chrUn_gl000225.fa,chrUn_gl000226.fa,chrUn_gl000227.fa,chrUn_gl000228.fa,chrUn_gl000229.fa,chrUn_gl000230.fa,chrUn_gl000231.fa,chrUn_gl000232.fa,chrUn_gl000233.fa,chrUn_gl000234.fa,chrUn_gl000235.fa,chrUn_gl000236.fa,chrUn_gl000237.fa,chrUn_gl000238.fa,chrUn_gl000239.fa,chrUn_gl000240.fa,chrUn_gl000241.fa,chrUn_gl000242.fa,chrUn_gl000243.fa,chrUn_gl000244.fa,chrUn_gl000245.fa,chrUn_gl000246.fa,chrUn_gl000247.fa,chrUn_gl000248.fa,chrUn_gl000249.fa,chrX.fa,chrY.fa human_hg19
```
You can also download the most used ones on the following link: **http://bowtie-bio.sourceforge.net/bowtie2/index.shtml**

Examples of lines used to launch the alignment procedure:
```bash
bowtie2 -x indices_genomes/sacCer3/sacCer3 -p6 --sam-no-hd --sam-no-sq --quiet --local --very-sensitive-local \
-S p1.sam /media/03b8b079-9d7a-4162-8201-6dd5d9923f62/2013/11_05_2013_Hi_Seq_MM/sequencage_nov2013/RSG6_L6/seq/BC76_CTGT.dat.end1.pcrfree

bowtie2 -x indices_genomes/sacCer3/sacCer3 -p6 --sam-no-hd --sam-no-sq --quiet --local --very-sensitive-local \
-S p2.sam /media/03b8b079-9d7a-4162-8201-6dd5d9923f62/2013/11_05_2013_Hi_Seq_MM/sequencage_nov2013/RSG6_L6/seq/BC76_CTGT.dat.end2.pcrfree
```
We could have also aligned the reads with an iterative alignment procedure like in [hiclib] (https://bitbucket.org/mirnylab/hiclib) ). This procedure start with a fixed lenght (i.e 20 bp), tries to align the read. If the read is correctly aligned, it is kept else the length is incremented (i.e by 5 bp) until a correct mapping can be found.


## Filtering of the data:
A removal of uncrosslinked events (uncuts, loops...) can be applied at this stage.  
This procedure is optional and might be necessary when you want to study the structure of chromatin at short scales like several kb. 

In this step, you need to assign a restriction fragment to every locus (chrm - position). 






## Session 2: Normalization of the data
We used the normalization procedure called SCN (presented in http://www.biomedcentral.com/1471-2164/13/436). 
This procedure assumes that every bin should be detected with the same strength. We divide each line by its sum then each column by its sum. We reiterate this loop several times. It has been shown that after several iterations, the matrice converged. 
Before the iterations, poor interacting bins must be discarded and considered as non detectable bins. 
To do that, a threshold is given as an argument to the SCN function so that every bin with a reads number below this value will be replaced by vectors of zeros. To determine the threshold you can plot the distribution in the number of reads by doing the histogram:


In the python code. 

It should be noticed that this procedure does not conserve the symetry propertie of the matrix. To 




### Session 3: Computation of a genomic distance law
This plot is very important and must be computed in the firts momemnt of the analysis. It reflects the polymer behavoir of the chromatin and allow to check if there is 3D signal. 
It consist in computing the mean number of reads in function of the genomic distance between the two loci. 


### Session 4: Computation of correlation matrices


### Session 5: Decomposition into eigen vectors 


### Session 6: Use of sparce formalism for contacts maps
Another way to mathematically represent the data, it the sparse formalism. It is very relevant for matrices in which most of the elements are zero which is ofter the case for human, mouse contacts maps. 








