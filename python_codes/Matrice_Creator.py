# -*- coding: utf-8 -*-
"""
To convert an alignment output file into a matrice object. 
author: Axel KournaK
"""
import numpy as np
import matplotlib
from pylab import *
import os
import sys

BIN= int(sys.argv[2]);    # Size (in bp) of the bin to build the contacts map 
mat={};   # dictionary object to keep in memory all the contacts
maxi={};  # dictionary object to keep in memory the lenghts of the chromosomes (in bins) 
i=0;

with open(sys.argv[1]) as f: # open the file for reading output alignment file
    for line in f: # iterate over each line
        i=i+1;
        if i % 1000000 == 0:
            print(str(i)+" lines parsed.");
        chr1, locus1, sens1,indice1, chr2, locus2, sens2,indice2 = line.split(); # split it by whitespace
        locus1=int(locus1);sens1=int(sens1);
        locus2=int(locus2);sens2=int(sens2); 
        bin1 = int(locus1 /  BIN);
        bin2 = int(locus2 /  BIN);
        key1=(chr1, bin1, chr2, bin2);
        key2=(chr2, bin2, chr1, bin1);
        
        if key1 in mat:
            mat[key1] += 1;
        else:
            mat[key1] = 1  
        if key2 in mat:
            mat[key2] += 1
        else:
            mat[key2] = 1;   
        if chr1 not in maxi:
            maxi[chr1] =0;  
        if chr2 not in maxi:
            maxi[chr2] =0;       
        if bin1 > maxi[chr1]:
            maxi[chr1] =  bin1;      
        if bin2 > maxi[chr2]:
            maxi[chr2] =  bin2;
 
# Check for the maximum of bins:  
list_chr = sort(maxi.keys() );   #  by default all the chromosomes of the organism
#list_chr = ["chr3"];  # or list of subset of chromosomes you want to display

N_BINS=0;
for chr in list_chr:
    N_BINS=N_BINS+ maxi[chr]+1;
    print(chr,maxi[chr]) ;
# 
print(N_BINS);
 
# Conversion of the Matrix in an numpy array object:
bin_mat1=0;
bin_mat2=0;
MATRICE=np.zeros( (N_BINS,N_BINS) );

for chr1 in list_chr :
    for bin1 in range(0,maxi[chr1]+1) :
        bin_mat2=0;
        for chr2 in list_chr  :  
            for bin2 in range(0,maxi[chr2]+1) :
                key1=(chr1, bin1, chr2, bin2);
                if key1 in mat:
                    mat[key1] = mat[key1]; 
                    MATRICE[bin_mat1,bin_mat2]= mat[key1];
                else:
                    mat[key1] = 0;
                    MATRICE[bin_mat1,bin_mat2]= 0;
                bin_mat2 +=1;
        bin_mat1 +=1;

# Plot of the matrice and savings:
imshow( MATRICE**0.2,interpolation="none");
colorbar();
savefig('MAT_RAW.png');
np.savetxt('MAT_RAW.txt',MATRICE)

