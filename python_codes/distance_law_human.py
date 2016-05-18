# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 14:58:56 2015
To compute a distance law plot from a Matrix.
@author: Axel
"""
#  DIST FUNCTION :

import numpy as np

def dist_law(A):
    n1 = A.shape[0];
    print "Size of the matrix entetered for the genomic distance law:"
    print n1;
    dist = np.zeros((n1, 1));
    for nw in range(0,n1) :   # scales
        somme = [];         
        for i in range(0,n1) :
                kp= i -nw;
                if kp >= 0 and  kp < n1 :       
                    somme.append(A[i,kp]);
        dist[nw] = np.mean(somme);                      
    return dist;


