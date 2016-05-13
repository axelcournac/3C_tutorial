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
We could have also aligned the reads with an iterative alignment procedure like in [hiclib] (https://bitbucket.org/mirnylab/hiclib) ). 


## Filtering of the data:
A removal of uncrosslinked events (uncuts, loops...) can be applied at this stage.  
This procedure is optional and might be necessary when you want to study the structure of chromatin at short scales like several kb. 





## Session 2: Normalization of the data
We used the normalization procedure called SCN (presented in http://www.biomedcentral.com/1471-2164/13/436). 
In the python code. 

It should be noted that this procedure does not conserv the symetry propertie of the matrix. 




### Session 3: Computation of a genomic distance law


### Session 4: Computation of correlation matrices, 


### Session 5: Decomposition into eigen vectors 








