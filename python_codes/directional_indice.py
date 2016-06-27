# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 16:46:59 2014
@author: axel KournaK
To calculate the directional preference for a given bin as in the Laub paper.
"""

from scipy import stats
import numpy as np 
import math

def directional(A, nw):
    n1 = A.shape[0]  
    print("Size of the matrix entetered for the directional index:")
    print(n1)
    signal1 = np.zeros((n1, 1));
    
    for i in range(0,n1) :
        vect_left = [];
        vect_right = [];
        
        for k in range(i-1,i-nw-1,-1) :
            kp =k; 
            if k < 0 :
                kp = n1 +k ;
            if A[i,kp] > 0 :
                vect_left.append(math.log(A[i,kp]));    
            else :
                vect_left.append(0);  
                    
                
        for k in range(i+1,i+nw+1) : 
            kp =k;
            if k >= n1 :
                kp = k - n1;
            if A[i,kp] > 0 :
                vect_right.append(math.log(A[i,kp]));    
            else :
                vect_right.append(0);  
                           
        if sum(vect_left) != 0 and sum(vect_right) != 0 :
            signal1[i] =  stats.ttest_rel(vect_right,vect_left)[0];
        else :
            signal1[i] =  0;
                           
    return signal1
    
#  for a bar graph :

#ind = np.arange(len(dir2))
#
#ind = np.arange(len(M)) 
#plt.bar(ind,M,color="red");
#
#dom22 =  [i/10 for i in dom22 ]
#
#ind = np.arange(len(dom22))
#plt.bar(ind, dom22)    
#show();
