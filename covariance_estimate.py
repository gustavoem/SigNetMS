# This modules defines a function that calculates the covariance matrix
# of variables given samples.

import numpy as np

def calc_covariance (sample):
    """ Compute an estimate of the covariance matrix of a normal
        distribution.

        Parameters
            sample: a list of list. Each element of the outmost list is
            a sample point of a multivariate normal random variable.
        
        Returns
            cov_matrix: an estimate of the covariance of the sample
            distribution.
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
    cov_matrix = cov_matrix / (len (sample) - 1)
    return cov_matrix
