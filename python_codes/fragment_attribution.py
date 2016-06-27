# -*- coding: utf-8 -*-
"""
This code assigns a restriction fragment to a locus (chrm-position). 
@author: Axel KournaK 
"""
import os
import sys
import glob
from pylab import *
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio import Restriction
from Bio.Seq import Seq
from Bio.Restriction import *

print("Arguments to enter:")
print("Path for fasta files of genome (ex: /run/media/axel/human_genome/), Restriction enzyme (ex: HindIII), Alignment file (ex: /data/output_alignment_idpt.dat)")

#path="/run/media/axel/RSG3/human_genome/";   #  path containing the fasta files of the chromosomes 
path=sys.argv[1]+"/";

list_files1 =  glob.glob(path+"*.fa");
list_files2 =  glob.glob(path+"*.fasta");
list_files3 =  glob.glob(path+"*.fas");
list_files4 =  glob.glob(path+"*.fna");
list_files=list_files1+list_files2+list_files3;
print("List of fasta files found:";
print(list_files; 

#enz=HindIII;   #   Restriction enzyme used in the experiment 
enz=sys.argv[2];
print("Enzyme used:"+str(enz))
print("Output file used:")
print(sys.argv[3];

rb = RestrictionBatch([enz]);    #  Restriction batch containing the restriction enzyme
restriction_table={};

for chr_file in list_files :
    handle = open(chr_file, "rU")
    for rec in SeqIO.parse(handle, "fasta"):
        print("Cutting the following chromosome:")
        chr_name=rec.id;
        print(chr_name)
        # Building of the restriction map:
        #map_restriction = enz.search(rec.seq);
        map_restriction = rb.search(rec.seq);
        enzyme = rb.get(enz);  #  conversion from string type to restriction type
        map_restriction = map_restriction[enzyme];
        
        start_positions= map_restriction;
        # Adding start position of the chromosome to the start positions of restriction fragments:
        start_positions = np.insert(start_positions, 0, 0);   
        
        end_positions= map_restriction; 
        # Adding end position of the chromosome to the end positions of restriction fragments:
        end_positions = np.insert(end_positions, len(end_positions), len(rec));
        restriction_table[chr_name] = vstack((start_positions,end_positions));
        print("Number of restriction sites for "+chr_name+" for the enzyme "+str(enzyme))
        print(len(map_restriction))

    handle.close();

print("Restriction maps of all chromosomes are now in memory!")

#  Function to assign the corresponding restriction fragment to a locus.
def find_frag(chrm,locus) :
    T=restriction_table[chrm];
    bottom = 0;            # indice of the first restriction fragment 
    top = len(T[1]) -1;   # indice of the last restriction fragment 
    while 1 : 
        index = int((top + bottom)/2);  # $index <- middle of the chrm
        if locus >= T[0,index] and locus < T[1,index]  :
            indice = index;
            break
        elif locus < T[0,index] :
            top = index - 1;
        elif locus >= T[1,index] :
            bottom = index + 1;    
    return indice

# Reading of the output file from the alignment 
fout = open(sys.argv[3]+".indices","w")   #  output file that will contain the indices 
#fout = open("/home/axel/Bureau/output_alignment_idpt.dat3","w");

i=0;
#with open("/home/axel/Bureau/output_alignment_idpt.dat") as f: # open the file for reading
with open(sys.argv[3]) as f: # open the file for reading
    for line in f: # iterate over each line
        i=i+1;
        if i % 1000000 == 0:
            print(str(i)+" lines have been attributed restriction fragments.";
        chr1, locus1, sens1, chr2, locus2, sens2 = line.split(); # split it by whitespace
        locus1=int(locus1);sens1=int(sens1);
        locus2=int(locus2);sens2=int(sens2);        
        indice1=find_frag(chr1,locus1);
        indice2=find_frag(chr2,locus2);
        fout.write(str(chr1)+"\t"+str(locus1)+"\t"+str(sens1)+"\t"+str(indice1)+"\t"+str(chr2)+"\t"+str(locus2)+"\t"+str(sens2)+"\t"+str(indice2)+"\n"); 

fout.close();
print("Output file done!")

