# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 11:08:24 2016
@author: Axel KournaK
To make histograms like in R. 
"""

import numpy as np
from pylab import *

def histo_r(V,number_of_bins):    
    hist, bins = np.histogram(V,bins= number_of_bins)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.show(); 
    return center,hist,width
