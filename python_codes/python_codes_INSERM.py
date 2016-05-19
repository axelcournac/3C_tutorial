import numpy as np
import matplotlib
from histo_r import *
import scn_human
import distance_law_human
import directional_indice
import matplotlib.gridspec as gridspec
import scipy as sc
from scipy import signal
from scipy.linalg import eig

# Plot of historgram of number of reads per bin:
histo_r(MATRICE.sum(axis=0),100);
savefig('chr3_histogram.png');

# Normalisation of the contacts map:
mn=scn_human.scn_func(MATRICE,1500);
imshow(mn**0.2,interpolation="none");
colorbar();
savefig('chr3_NORMALISED.png');
np.savetxt('chr3_NORMALISED.txt',MATRICE);


# Computation and plot of genomic distance law:
d=distance_law_human.dist_law(mn);
plt.plot(d,label="Chr 3",linewidth=3.0);
plt.loglog();
plt.xlabel("Genomic distance in bins of 200kb");
plt.ylabel("Contacts Frequency");
plt.legend();
savefig('chr3_distance_law.png');

# Computation of the correlation matrice
n1=mn.shape[0];
# Computation of genomic distance law matrice:
MAT_DIST =  np.zeros((n1, n1));
for i in range(0,n1) :
    for j in range(0,n1) :
        MAT_DIST[i,j] =  d[abs(j-i)] ;
        
#imshow(MAT_DIST**0.2);
MAT2=mn/MAT_DIST;
imshow( MAT2**0.2,interpolation="none");

# Computation of the correlation matrice:
mc=np.corrcoef(MAT2);
mc[np.isnan(mc)] = 0.0;
imshow(  (mc-mc.min())**0.2,interpolation='none');
savefig('chr3_CORRELATION.png');

mc2=np.corrcoef(mc);
mc2[np.isnan(mc2)] = 0.0;

imshow( mc2,interpolation='none');
colorbar();
savefig('chr3_CORRELATION2.png');

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
