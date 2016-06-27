# lines used to use Jupyter

import numpy as np
import os
import matplotlib
%matplotlib inline
import mpld3
mpld3.enable_notebook()

pylab.rcParams['figure.figsize'] = (10, 10)
m=loadtxt("/data/tpdata/hic/MAT_RAW.txt")
imshow(m**0.2,interpolation="none")
