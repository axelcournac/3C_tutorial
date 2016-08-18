# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from pylab import *  
    
# Normalisation function
# A is the matrix to normalise  
# the treshold to remove bins for the normalization (these bins are replaced by zeros vectors)
     
def scn_func(A,threshold):    
    n1 = A.shape[0];
    n_iterations=10;
    keep = np.zeros((n1, 1));
    
    for i in range(0,n1) :
        if np.sum(A[i,]) > threshold:
            keep[i] = 1
        else :
            keep[i] = 0
    
    indices1=np.where(keep >0 )
    indices2=np.where(keep <=0 )
    
    for n in range(0,n_iterations) :
        print(n);
        for i in range(0,n1) :
            A[indices1[0],i]=A[indices1[0],i]/ np.sum(A[indices1[0],i])
            A[indices2[0],i]=0   
        A[np.isnan(A)] = 0.0 
        
        for i in range(0,n1) :    
            A[i,indices1[0]]=A[i,indices1[0]]/ np.sum(A[i,indices1[0]])
            A[i,indices2[0]]=0  
        A[np.isnan(A)] = 0.0    
        
    return A
