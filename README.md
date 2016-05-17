# Tutorial for 3C/Hi-C data processing 


This repository contains several codes and simple command lines for the processing, vizualisation and primary analysis of 3C/Hi-C datasets.
These codes were used during the INSERM workshop "Capturing chromosone conformation: toward a 3D view of genome regulation", May 9-13, Paris at Pasteur Institut.

For queries or help getting these running, you can contact me on mail or open an issue at the github repository.

### Dependencies

Scripts will run on OS X and other Unix-based systems. External dependencies should be installed somewhere on your `$PATH`.
It basically requires to have Python installed on your machine which is commonly installed on Unix-based systems. 
For windows, you can have a look to https://www.python.org/downloads/windows/. Then, a few specific python modules are necessary for diverses operations on arrays and vizualisation. 

#### Python
* Numpy
* Matplotlib
* Scipy
* Biopython

#### External programs

* `SRA tool` / [SRA](http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software)
* `Bowtie2 ` / [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

## Raw data extraction and alignment
#### Data extraction
Raw data are deposited on Short Read Archive server at the following address **http://www.ncbi.nlm.nih.gov/sra**.
In this tutorial, we take as example one run of the replicat 2 of IMR90 from [Dixon et al. Nature 2012] (http://www.nature.com/nature/journal/v485/n7398/full/nature11082.html). Raw data can be retreived at  http://www.ncbi.nlm.nih.gov/sra/SRX212173 (the identification number is SRR639031). 

We used an SRA executable called fastq-dump from SRA to extract and split both mates of a library (to use it, you can go with your terminal to the directory containg the executables files by using the bash command cd).Then the program can be used like this:  /fastq-dump library_identification --split-3 -O /path_to_a_directory

```bash
./fastq-dump SRR639031 --split-3 -O /run/media/axel/RSG3/IMR90_data/
```


#### Alignment

Before aligning the reads on a reference genome, you need to index it, the line to enter is as follows: bowtie2-build [list of fasta files of your genome separated by comma] [name of index you want to give].

```bash
bowtie2-build chr10.fa,chr11.fa,chr11_gl000202_random.fa,chr12.fa,chr13.fa,chr14.fa,chr15.fa,chr16.fa,chr17_ctg5_hap1.fa,chr17.fa,chr17_gl000203_random.fa,chr17_gl000204_random.fa,chr17_gl000205_random.fa,chr17_gl000206_random.fa,chr18.fa,chr18_gl000207_random.fa,chr19.fa,chr19_gl000208_random.fa,chr19_gl000209_random.fa,chr1.fa,chr1_gl000191_random.fa,chr1_gl000192_random.fa,chr20.fa,chr21.fa,chr21_gl000210_random.fa,chr22.fa,chr2.fa,chr3.fa,chr4_ctg9_hap1.fa,chr4.fa,chr4_gl000193_random.fa,chr4_gl000194_random.fa,chr5.fa,chr6_apd_hap1.fa,chr6_cox_hap2.fa,chr6_dbb_hap3.fa,chr6.fa,chr6_mann_hap4.fa,chr6_mcf_hap5.fa,chr6_qbl_hap6.fa,chr6_ssto_hap7.fa,chr7.fa,chr7_gl000195_random.fa,chr8.fa,chr8_gl000196_random.fa,chr8_gl000197_random.fa,chr9.fa,chr9_gl000198_random.fa,chr9_gl000199_random.fa,chr9_gl000200_random.fa,chr9_gl000201_random.fa,chrM.fa,chrUn_gl000211.fa,chrUn_gl000212.fa,chrUn_gl000213.fa,chrUn_gl000214.fa,chrUn_gl000215.fa,chrUn_gl000216.fa,chrUn_gl000217.fa,chrUn_gl000218.fa,chrUn_gl000219.fa,chrUn_gl000220.fa,chrUn_gl000221.fa,chrUn_gl000222.fa,chrUn_gl000223.fa,chrUn_gl000224.fa,chrUn_gl000225.fa,chrUn_gl000226.fa,chrUn_gl000227.fa,chrUn_gl000228.fa,chrUn_gl000229.fa,chrUn_gl000230.fa,chrUn_gl000231.fa,chrUn_gl000232.fa,chrUn_gl000233.fa,chrUn_gl000234.fa,chrUn_gl000235.fa,chrUn_gl000236.fa,chrUn_gl000237.fa,chrUn_gl000238.fa,chrUn_gl000239.fa,chrUn_gl000240.fa,chrUn_gl000241.fa,chrUn_gl000242.fa,chrUn_gl000243.fa,chrUn_gl000244.fa,chrUn_gl000245.fa,chrUn_gl000246.fa,chrUn_gl000247.fa,chrUn_gl000248.fa,chrUn_gl000249.fa,chrX.fa,chrY.fa hg19
```
You can also download the most used ones on the following link: **http://bowtie-bio.sourceforge.net/bowtie2/index.shtml**

We align the raw reads with the software bowtie2 with several options. Information can be found in the following reference manual (**http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml**).

Examples of lines used to launch the alignment procedure:
```bash
bowtie2 -x /run/media/axel/Pasteur_Koszul/human_index_downloaded/hg19 -p6 --sam-no-hd --sam-no-sq --quiet --local --very-sensitive-local -S p1.sam SRR639031_1.fastq 

bowtie2 -x /run/media/axel/Pasteur_Koszul/human_index_downloaded/hg19 -p6 --sam-no-hd --sam-no-sq --quiet --local --very-sensitive-local -S p2.sam SRR639031_2.fastq 
```
We could have also aligned the reads with an iterative alignment procedure like in [hiclib] (https://bitbucket.org/mirnylab/hiclib) ). In this procedure, each read starts with a fixed lenght (i.e 20 bp), tries to align the read. If the read is correctly aligned, it is kept else the length is incremented (i.e by 5 bp) until a correct mapping can be found.
It is important to align each mate idependently and then repair them  (Bowtie expects a certain distribution of distances between mates so the "pairs mode" of Bowtie is not recommanded for Hi-C data). 
Here, some lines that can be used to do this task:

```bash
#  Keeping only the columns of the sam file that contain necessary information:
awk '{print $1,$3,$4,$2,$5;}' p1.sam > p1.sam.0
awk '{print $1,$3,$4,$2,$5;}' p2.sam > p2.sam.0

# Sort according to the read identification to have both mates in the same order
sort -d -k1 p1.sam.0 > p1.sam.0.sorted
sort -d -k1 p2.sam.0 > p2.sam.0.sorted

# Pairing of both mates in a single file
paste p1.sam.0.sorted p2.sam.0.sorted > p1_p2_merged

# Removal of intermediar files
rm p1.sam.0.sorted
rm p1.sam.0.sorted

# Filtering of paires of reads that both have a Mapping Quality above 30
awk '{if($1==$6 && $5>= 30 && $10 >= 30) print $2,$3,$4,$7,$8,$9}'  p1_p2_merged  > output_alignment_idpt.dat

# Removal of intermediar file
rm p1_p2_merged
```
These lines can be combined in a same bash script and be executed by launching the script:
```bash
bash script_bash_INSERM.bh
```

#### Filtering of the data
A removal of uncrosslinked events (uncuts, loops...) can be applied at this stage.  
This procedure is optional and might be necessary when you want to study the structure of chromatin at short scales like several kb. 

In this step, you need to assign a restriction fragment to every locus (chrm - position). This can be done like in the python code fragment_attribution.py. 

## Building of the contacts map

From the contacts, 


## Normalization of the data
We used the normalization procedure called SCN (for Sequential Component Normalization presented in http://www.biomedcentral.com/1471-2164/13/436). 
This procedure assumes that every bin should be detected with the same strength. We divide each line by its sum then each column by its sum. We reiterate this loop several times. It has been shown that after several iterations, the matrice converged. 
Before the iterations, poor interacting bins must be discarded and considered as non detectable bins. 
To do that, a threshold is given as an argument to the SCN function so that every bin with a reads number below this value will be replaced by vectors of zeros. To determine the threshold you can plot the distribution in the number of reads by doing the histogram:

```python
from histo_r import *
histo_r(m.sum(axis=0),100)
```
The function histo_r(V,N) makes a histogram of the vector V in N bins and plots the results as a bar plot (it is an equivalent of the R function hits).   

In the python code. 

It should be noticed that this procedure does not conserve the symetry propertie of the matrix. To recover the symetry property you can use a simple line in python:

```python
mn=(mn+mn.T)/2;
```
mn.T is the transposed matrice of mn. 



## Computation of genomic distance law
This plot is important and must be computed at the very first steps of data processing. It reflects the polymer behaviour of the chromatin and thus allows to check the presence or absence of 3D signal. 
It consists in computing the mean number of reads in function of the genomic distance separating the two loci. 
The computation thus consists in scanning every diagonal of the matrice and taking the average of this sets of elements.


## Computation of correlation matrices
The correlation matrice is very often computed in Hi-C data analysis. Even if not present in final publication, it can give a more visible image of the structures that can be detected especially domains structures (squares in the contacts maps). It consists in looking for correlation in the contacts patterns between each line and column. 
It is simply computed by taking the Pearson Coefficient (or another correlation coefficient) between each line and each column of a matrice. 


## Directional Index tool to detect TADs
This tool is commonly used in Hi-C data analysis. It looks for change in the directionality between "left vector" and "right vector" at a certain loci in the genome. A change could come from the presence of a border between two different compartments in the genome.
It consists in doing a paired T test between "left vector" and "right vector" on each bin along the genome. The size of the "left vector" and "right vector" is put as a parameter and allows to look for domains structures at a specific scale. 


## Decomposition into eigen vectors 
This geometrical transformation allows to decompose the matrice into eigen values and vectors. It is a way to simplify the data or at least to decrease the dimentionality of the mathematical object. 
It has been shown that the first eigen vector corresponds to the 2 compartments partition of the genome. 
It should be kept in mind that this is a first approximation of the 3D structure of a genome and other or sub-compartments can be detected using higher resolution and/or using other tools of compartment detection.


## Use of sparse formalism
An alternative and interesting way to mathematically represent the data is the sparse formalism. It is very relevant for matrices in which most of the elements are zero which is ofter the case for human, mouse contacts maps. 
![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/sparse.png)

In python, functions implemented for the use of sparse formalism can be found in the module scipy.

