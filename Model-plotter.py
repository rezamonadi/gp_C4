def build_correlation_matrix(M):
        '''
        Covert covariance matrix to correlation matrix

        Parameters:
        ----
        M (N pixels, k) : K = M M' , K is a covariance matrix, M is its matrix decomposition   

        Return:
        ----
        C : correlations matrix with its diag C = I
        '''
        # build covariance matrix
        K = np.matmul( M, M.T )

        # query diag elements
        d = np.sqrt(np.diag( K ))[:, np.newaxis]
        
        M_div_d = M / d

        C = np.matmul( M_div_d, M_div_d.T)

        return C

from scipy.io import loadmat                                                                            
import matplotlib.pyplot as plt
import numpy as np
M = loadmat("M.mat")                                                                                   
# M = loadmat("MM-70%-1350-1570.mat")                                                                                   
C = build_correlation_matrix(M["M"])    
min_lambda = 1216
max_lambda = 1600
plt.imshow(C, origin='low')
plt.show()
# locs, labels= plt.xticks()
# locs
# new_labels = np.linspace(min_lambda, max_lambda,len(locs)-1)
# new_labels = np.int16(new_labels)
# print(new_labels)
# x_ticks=[]
# plt.xticks(locs, new_labels)
# plt.savefig('MM-0.5-1350-1600.png')
plt.savefig('MM-1216-1600.png')