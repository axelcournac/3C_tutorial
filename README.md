# Tutorial for 3C/Hi-C data processing 


This document presents several codes and simple command lines for the processing, vizualisation and primary analysis of 3C/Hi-C datasets.
These codes were used during the INSERM workshop "Capturing chromosone conformation: toward a 3D view of genome regulation", May 9-13, Paris at Pasteur Institut.

For queries or help getting these running, you can send email or open an issue at the github repository.

### Table of contents

* [Dependencies](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#dependencies)
* [Raw data extraction and alignment](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#raw-data-extraction-and-alignment)
* [Building of the contacts map](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#building-of-the-contacts-map)
* [Normalization of the data](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#normalization-of-the-data)
* [Computation of genomic distance law](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#computation-of-genomic-distance-law)
* [Computation of correlation matrices](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#computation-of-correlation-matrices)
* [Directional Index tool to detect TADs](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#directional-index-tool-to-detect-tads)
* [Decomposition into eigen vectors](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#decomposition-into-eigen-vectors)
* [Use of sparse formalism](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#use-of-sparse-formalism)
* [Miscellaneous](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#miscellaneous)

### Dependencies

Scripts will run on OS X and other Unix-based systems. It basically requires to have Python installed on your machine which is commonly installed on Unix-based systems. 
For windows, you can have a look to https://www.python.org/downloads/windows/. Then, a few python modules are necessary for diverses operations on arrays and vizualisation. 

#### Python (>=2.7)
* Numpy
* Matplotlib (>=1.0)
* Scipy
* Biopython

#### External programs

* `SRA tool` / [SRA](http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software)
* `Bowtie2 ` / [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

## Raw data extraction and alignment
#### Data extraction
As example, we will use raw data that are deposited on Short Read Archive server at the following address **http://www.ncbi.nlm.nih.gov/sra**.
We take as example one run of the replicat 2 of IMR90 from [Dixon et al. Nature 2012] (http://www.nature.com/nature/journal/v485/n7398/full/nature11082.html). Raw data can be retreived at  http://www.ncbi.nlm.nih.gov/sra/SRX212173 (the identification number is SRR639031). 

We used an SRA executable called fastq-dump from SRA to extract and split both mates of a library (to use it, you can go with your terminal to the directory containg the executables files by using the bash command cd).Then the program can be used like this:  /fastq-dump library_identification --split-3 -O /path_to_a_directory

```bash
./fastq-dump SRR639031 --split-3 -O /run/media/axel/RSG3/IMR90_data/
```


#### Alignment

For the alignement step, we will use the sofware Bowtie2. Before aligning the reads on a reference genome, it is necessary to index it, the line to enter is as follows: bowtie2-build [list of fasta files of your genome separated by comma] [name of index you want to give].

```bash
bowtie2-build chr10.fa,chr11.fa,chr11_gl000202_random.fa,chr12.fa,chr13.fa,chr14.fa,chr15.fa,chr16.fa,chr17_ctg5_hap1.fa,chr17.fa,chr17_gl000203_random.fa,chr17_gl000204_random.fa,chr17_gl000205_random.fa,chr17_gl000206_random.fa,chr18.fa,chr18_gl000207_random.fa,chr19.fa,chr19_gl000208_random.fa,chr19_gl000209_random.fa,chr1.fa,chr1_gl000191_random.fa,chr1_gl000192_random.fa,chr20.fa,chr21.fa,chr21_gl000210_random.fa,chr22.fa,chr2.fa,chr3.fa,chr4_ctg9_hap1.fa,chr4.fa,chr4_gl000193_random.fa,chr4_gl000194_random.fa,chr5.fa,chr6_apd_hap1.fa,chr6_cox_hap2.fa,chr6_dbb_hap3.fa,chr6.fa,chr6_mann_hap4.fa,chr6_mcf_hap5.fa,chr6_qbl_hap6.fa,chr6_ssto_hap7.fa,chr7.fa,chr7_gl000195_random.fa,chr8.fa,chr8_gl000196_random.fa,chr8_gl000197_random.fa,chr9.fa,chr9_gl000198_random.fa,chr9_gl000199_random.fa,chr9_gl000200_random.fa,chr9_gl000201_random.fa,chrM.fa,chrUn_gl000211.fa,chrUn_gl000212.fa,chrUn_gl000213.fa,chrUn_gl000214.fa,chrUn_gl000215.fa,chrUn_gl000216.fa,chrUn_gl000217.fa,chrUn_gl000218.fa,chrUn_gl000219.fa,chrUn_gl000220.fa,chrUn_gl000221.fa,chrUn_gl000222.fa,chrUn_gl000223.fa,chrUn_gl000224.fa,chrUn_gl000225.fa,chrUn_gl000226.fa,chrUn_gl000227.fa,chrUn_gl000228.fa,chrUn_gl000229.fa,chrUn_gl000230.fa,chrUn_gl000231.fa,chrUn_gl000232.fa,chrUn_gl000233.fa,chrUn_gl000234.fa,chrUn_gl000235.fa,chrUn_gl000236.fa,chrUn_gl000237.fa,chrUn_gl000238.fa,chrUn_gl000239.fa,chrUn_gl000240.fa,chrUn_gl000241.fa,chrUn_gl000242.fa,chrUn_gl000243.fa,chrUn_gl000244.fa,chrUn_gl000245.fa,chrUn_gl000246.fa,chrUn_gl000247.fa,chrUn_gl000248.fa,chrUn_gl000249.fa,chrX.fa,chrY.fa hg19
```
You can also download the most used indexes on the following link: **http://bowtie-bio.sourceforge.net/bowtie2/index.shtml**

We align the raw reads with the software bowtie2 with several options. Information can be found in the following reference manual (**http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml**).

Examples of lines used to do the alignment procedure:
```bash
bowtie2 -x /run/media/axel/human_index_downloaded/hg19 -p6 --sam-no-hd --sam-no-sq --quiet --local --very-sensitive-local -S p1.sam SRR639031_1.fastq 

bowtie2 -x /run/media/axel/human_index_downloaded/hg19 -p6 --sam-no-hd --sam-no-sq --quiet --local --very-sensitive-local -S p2.sam SRR639031_2.fastq
```
Remark: We could have also aligned the reads with an iterative alignment procedure like in [hiclib] (https://bitbucket.org/mirnylab/hiclib). In this procedure, each read starts with a fixed lenght (i.e 20 bp), the aligner tries to find a correct (unambiguous) location. If the read can be correctly aligned, it is kept else the length is incremented (i.e by 5 bp) until a correct mapping can be found.
It is important to align each mate idependently and then repair them  (Bowtie expects a certain distribution of distances between mates so the "pairs mode" of Bowtie is not recommanded for Hi-C data). 
Here, some lines that can be used to do this latter task:

```bash
#  Keeping only the columns of the sam file that contain necessary information:
awk '{print $1,$3,$4,$2,$5;}' p1.sam > p1.sam.0
awk '{print $1,$3,$4,$2,$5;}' p2.sam > p2.sam.0

# Sort according to the read identification to have both mates in the same order
# if sort does not have -V option try -d
sort -V -k1 p1.sam.0 > p1.sam.0.sorted
sort -V -k1 p2.sam.0 > p2.sam.0.sorted

# Pairing of both mates in a single file
paste p1.sam.0.sorted p2.sam.0.sorted > p1_p2_merged

# Removal of intermediar files
rm p1.sam.0.sorted
rm p2.sam.0.sorted

# Filtering of paires of reads that both have a Mapping Quality above 30
awk '{if($1 eq $6 && $5>= 30 && $10 >= 30) print $2,$3,$4,$7,$8,$9}'  p1_p2_merged  > output_alignment_idpt.dat

# Removal of intermediar file
rm p1_p2_merged
```
These lines can be combined in a same bash script and be executed by launching the script [`script_bash_INSERM.bh`](bash_lines/script_bash_INSERM.bh):
```bash
bash script_bash_INSERM.bh
```

At this stage, you should have a file containing these information:
```
chrX 104115113 16 chr5 169107262 0 
chr15 64627253 16 chr15 64627696 0 
chr1 9155504 16 chr10 77551370 16 
chr20 36390507 0 chr15 48945063 0 
chrX 134656481 0 chrX 134656772 16 
chr4 127660544 16 chr1 7470584 0 
chr5 137488897 0 chr9 27468101 0 
chr3 178179845 0 chr3 171257372 16
```

The 3rd and 6th fields correspond to the directions of the reads: 0 means in the same direction of the reference genome, 16 means in the opposite direction.
You can also assign a restriction fragment to each locus using the python code [`fragment_attribution.py`](python_codes/fragment_attribution.py):

The command to enter to use it is: 

python fragment_attribution.py [path of fasta files of the genome] [Restriction enzyme]  [Alignment file]

Ex: 
```bash
python fragment_attribution.py /run/media/axel/RSG3/human_genome/ HindIII /run/media/axel/RSG3/IMR90_data/output_alignment_idpt.dat
```
The indices start at 0 for every chromosome, you should have a file now containing the indices in 4th and 8th columns:
```
chrX	104115113	16	29877	chr5	169107262	0	52159
chr15	64627253	16	12632	chr15	64627696	0	12632
chr1	9155504	16	1806	chr10	77551370	16	21377
chr20	36390507	0	8974	chr15	48945063	0	7863
chrX	134656481	0	39113	chrX	134656772	16	39113
chr4	127660544	16	38703	chr1	7470584	0	1411
chr5	137488897	0	42558	chr9	27468101	0	8715
chr3	178179845	0	53834	chr3	171257372	16	51648
```

#### Filtering of the data
A check for the proportion of uncrosslinked events (uncuts, loops...) can be carried out at this stage.  
With these information, one can remove uninformative events more cautiously then just removing all events below 10 kb. This filtering is optional and might be necessary when you want to study the structure of chromatin at very short scales like several kb. 
We used the home-made python code [`library_events.py`](python_codes/library_events.py):
```bash
python library_events.py /run/media/axel/RSG3/IMR90_data/output_alignment_idpt.dat.indices
```

You can have a look to the behaviors of the different pairs of reads taking into account the directions of the reads with this kind of plots:

![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/behavior_events2.png)

At this step, close the window, the program will ask you the two thresholds for the uncuts and loops events that you can determine from the previous plot.
Then, the code counts the different types of events, it gives information about the quality of your library preparation.
![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/PieChart_events3.png)


## Building of the contacts map

From the alignment file, you can build a binned contacts map. We used the home-made python code [`Matrice_Creator.py`](python_codes/Matrice_Creator.py): that converts the alignment file information into a dictionary structure and then into an array of the chosen chromosome(s). In this code, there are 3 parameters you can modify: 
the size of the bin, the path of your output alignment file and the list of the chromosome(s) you want to display.

For the human chromosome 3, you shoud have a picture looking like this:

![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/chr3_RAW.png)

## Normalization of the data
We used the normalization procedure called SCN (for Sequential Component Normalization presented in http://www.biomedcentral.com/1471-2164/13/436). 
This procedure assumes that every bin should be detected with the same strength. We divide each line of the previous matrice by its sum then each column by its sum. We reiterate this loop several times. It has been shown that after several iterations, the matrice has converged. 
Before the iterations, poor interacting bins must be discarded and considered as non detectable bins (due to experimental biases, unmapple sequences etc). 
To do that, a threshold is given as an argument to the SCN function so that every bin with a reads number below this value will be replaced by vectors of zeros. To determine this threshold, you can plot the distribution in the number of reads by doing the histogram in a python terminal: 

```python
from pylab import *
from histo_r import *

MATRICE=loadtxt("chr3_RAW.txt");
histo_r(MATRICE.sum(axis=0),100);
```
The function histo_r(V,N) makes a histogram of the vector V in N bins and plots the result as a bar plot (it is an equivalent of the R function hist() ).  

![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/histogram_chr3.png)

To normalise the raw contact map, we used the function scn_func coded in the module [scn_human](https://github.com/axelcournac/3C_tutorial/blob/master/python_codes/scn_human.py). 
```python
import scn_human

mn=scn_human.scn_func(MATRICE,1500);
imshow(mn**0.2,interpolation="none");
colorbar();
savefig('chr3_NORMALISED.png');
np.savetxt('chr3_NORMALISED.txt',MATRICE);
```
The normalised contacts map should look like this:

![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/chr3_NORMALISED.png)

The matrice elements can be viewed as relative frequencies of contacts between different loci (the scores go now from 0 to 1).
It should be noticed that this procedure does not conserve the symetry property of the matrix. To recover the symetry property, you can use a simple line in python:
```python
mn=(mn+mn.T)/2;
```
where mn.T is the transposed matrice of mn. 

## Computation of genomic distance law
This plot is important and must be computed at the very first steps of data processing. It reflects the polymer behaviour of the chromatin and thus allows to check the presence or absence of 3D signal. 
It consists in computing the mean number of reads in function of the genomic distance separating the two loci. 
The computation thus consists in scanning every diagonal of the matrice and taking the average of this sets of elements.
We coded this computation in the function [distance_law_human](https://github.com/axelcournac/3C_tutorial/blob/master/python_codes/distance_law_human.py). 

```python
from pylab import *
import distance_law_human

# mn: normalised matrice
d=distance_law_human.dist_law(mn)
plt.plot(d,label="Chr 3",linewidth=3.0);
plt.loglog();
plt.xlabel("Genomic distance in bins of 100kb");
plt.ylabel("Contacts Frequency");
plt.legend();
savefig('chr3_distance_law.png');
```
Example of plot obtained for the chromosome 3:

![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/chr3_distance_law.png)



## Computation of correlation matrices
The correlation matrice is often computed in Hi-C data analysis. Even if not always present in final publication, it can render  more visible structures that are present in the data especially domains structures (squares in the contacts maps). It consists in looking for correlation in the shared neighbors between each line and column. 
It is simply computed by taking the Pearson Coefficient (or another correlation coefficient) between each line and each column of a matrice coded in the function corrcoef from Numpy.

```python
n1=mn.shape[0];
# Computation of genomic distance law matrice:
MAT_DIST =  np.zeros((n1, n1));
for i in range(0,n1) :
    for j in range(0,n1) :
        MAT_DIST[i,j] =  d[abs(j-i)] ;
        
#imshow(MAT_DIST**0.2);
MAT2=mn/MAT_DIST;
imshow( MAT2**0.2);

# Computation of the correlation matrice:
mc=np.corrcoef(MAT2);
mc[np.isnan(mc)] = 0.0;

mc2=np.corrcoef(mc);
mc2[np.isnan(mc2)] = 0.0;

imshow( mc2,interpolation='none');
colorbar();
savefig('chr3_CORRELATION2.png');
```

Example of plot for the human chromosome 3 with 100 kb bins and 2 iterations of correlation:
![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/chr3_CORRELATION2.png)


## Directional Index tool to detect TADs
This tool is commonly used in Hi-C data analysis. It looks for change in the directionality between "left vector" and "right vector" at a certain loci in the genome. A change could come from the presence of a border between two different compartments in the genome.
It consists in doing a paired T test between "left vector" and "right vector" on each bin along the genome. The size of the "left vector" and "right vector" is put as a parameter and allows to look for domains structures at a specific scale. 
We coded this computation in the function [directional_indice](https://github.com/axelcournac/3C_tutorial/blob/master/python_codes/directional_indice.py). 

```python
# Computation and plot of Directional Index:
M = np.corrcoef(mn);
n1 = M.shape[0];
scale=20;
DI = directional_indice.directional(M,scale).T;

# PLOT of contacts Map and Directional Index:
gs = gridspec.GridSpec(2, 1, height_ratios=[5,1] );  #  to have several subplots on the same graph
ax0 = plt.subplot(gs[0]);
ax0.imshow(mn**0.2,interpolation="none")
ax0.set_title("Contacts Map");

ax2 = plt.subplot(gs[1], sharex=ax0);
b1=0;
b2=len(DI.T);
borders = list(range(b1,b2));
borders2 = DI.T 
borders2 = array(borders2[:,0]);

ax2.set_xlim([b1, b2]);
ax2.set_ylim([-2.0, 2.0]);
ax2.fill_between( borders, 0,borders2, color='red' );
ax2.fill_between( borders, 0,borders2,borders2<0 ,color='green' );
ax2.set_ylabel("DI scale = 2Mb)");
savefig('chr3_DI_bin100kb.png');
```


Example of plot for the human chromosome 3 with 100 kb bins and DI carried out at the scale of 2Mb:
![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/chr3_DI_bin100kb.png)



## Decomposition into eigen vectors 
This geometrical transformation allows to decompose the matrice into eigen values and vectors. It is a way to simplify the data or at least to decrease the dimentionality of the mathematical object. 
It has been shown that the first eigen vector corresponds to the 2 compartments partition signal of the genome. 
It should be kept in mind that this is a first approximation of the 3D structure of a genome and other or sub-compartments can be detected using higher resolution and/or using other tools of compartment detection.

```python
# Computation of eigen vectors:
(V,D) = sc.linalg.eig(mc);
plot(D[:,0]);

# Plot of contacts Map and Eigen vector:
gs = gridspec.GridSpec(2, 1, height_ratios=[5,1] );
ax0 = plt.subplot(gs[0]);
ax0.imshow(mn**0.2,interpolation="none")
ax0.set_title("Contacts Map");

ax2 = plt.subplot(gs[1], sharex=ax0 );
b1=0;
b2=len(D[:,0]);
borders = list(range(b1,b2));
borders2 = D[:,0] 
borders2 = array(borders2[:,0]);

ax2.set_xlim([b1, b2]);
#ax2.set_ylim([-2.0, 2.0]);
ax2.fill_between( borders, 0,borders2, color='darkblue');
ax2.fill_between( borders, 0,borders2,borders2<0 ,color='darkblue' );
ax2.set_xlabel("Position along the genome (in bins of 200kb)");
ax2.set_ylabel("Eigen vector");
savefig('chr3_Eigen_bin100kb.png');
```


Example of plot for the human chromosome 3 with 100 kb bins (zoom at the end of the chromosome):
![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/chr3_Eigen_bin100kb.png)

## Use of sparse formalism
An alternative and interesting way to mathematically represent the data is the sparse formalism. It is very relevant for matrices in which most of the elements are zero which is ofter the case for human, mouse contacts maps. 
![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/sparse.png)

In python, functions implemented for the use of sparse formalism can be found in the module scipy.


## Miscellaneous
From the contacts map, you can do a kind of 4C plot of a position of interest by simply select the corresponding line in the matrix:

```python
plot(mn[ int(37000000/BIN),: ] );
xlabel("Postion along the chromosome in bins of 100kb")
ylabel("Normalised Score");
savefig('4Cplot_chr3.png');
```
Example of 4C plot from the normalised contacts map for the position 37000000 bp on the chromosome 3:
![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/4Cplot_chr3.png)







