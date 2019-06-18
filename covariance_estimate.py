# This modules defines a function that calculates the covariance matrix
# of variables given samples.

import numpy as np

def calc_covariance_diagonal (sample):
    """ Given a list of samples (lists of values) of variables, compute
        an estimation of a matrix that has only the diagonal of the
        covariance matrix and zeros on other positions. """
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
    return cov_matrix.diagonal () * np.eye (n)

     
