# This modules defines a function that calculates the covariance matrix
# of variables given samples.

import numpy as np

def calc_covariance (sample):
    """ Compute diagonal of covariance of a sample.

        Given a list of samples (lists of values) of variables, compute
        an estimation of the covariance matrix and return a matrix that 
        only has its diagonal.
        
        Parameters
            sample: a list of list. The inner list is a sample point of
            a multivariate random variable.
        
        Returns
            cov_matrix: a diagonal matrix M such that (M)i,i is the 
            variance of the i-th component of the random variable. This
            matrix is a numpy array object.
        """
    n = len (sample[0])
    theta_sum = np.zeros (n)
    for theta in sample:
        theta_sum += theta
    theta_mean = theta_sum / len (sample)
       
    cov_matrix = np.zeros ([n, n])
    for j in range (len (sample)):
        thetaj = sample[j]
        v = thetaj - theta_mean
        v = v.reshape ([n, 1])
        vT = np.array (v.transpose ())
        cov_matrix += v * vT
    cov_matrix = cov_matrix / len (sample)
    return cov_matrix
