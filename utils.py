import matplotlib
matplotlib.use('Agg')

from lxml import etree
import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import errno

def safe_log (x):
    """ Calculates the log of x safely.
        
        Parameters
            x: a floating point.

        Returns
            np.log (x), if x is greater than the smallest float 
                representable, and -inf otherwise.
    """
    if x < 0:
        raise ValueError ("x needs to be nonnegative.")

    min_positive = sys.float_info.min
    if x < min_positive:
        return float ("-inf")
    else:
        return np.log (x)


def safe_power (a, b):
    """ Calculates a ** b safely.
    
        Parameters
            a: a number
            b: also a number

        Returns
            a ** b if that number is smallest than the biggest possible
                float and inf otherwise.
    """
    try:
        x = a ** b
    except OverflowError:
        x = float ("inf")
    return x


def safe_pow_exp_ratio (a, b, t): 
    """ Calculates the ratio pow (exp (a) / exp (b), t). 
    
    Parameters
        a: a float.
        b: a float.
        t: a float.

    Returns
        +inf if a is not -inf and b is -inf
        -inf if a is -inf and b is not -inf
        1 if a is -inf and b is -inf (this is not true, I know, but it 
            is useful for calculating jumping probability rations in 
            MCMC).
        (exp (a) / exp (b))^t otherwise.
    """
    # if a is not -inf and b is -inf: I want the ratio to be +inf
    # if a is -inf and b is not -inf: I want the ratio to be -inf
    # if a is -inf and b is -inf: I want the ratio to be 1
    if not a > float ("-inf"):
        at = a
    else:
        at = a * t

    if not b > float ("-inf"):
        bt = b
    else:
        bt = b * t
    
    return safe_exp_ratio (at, bt)


def safe_exp_ratio (a, b):
    """ Calculates the ratio exp (a) / exp (b). 
    
    Parameters
        a: a float.
        b: a float.

    Returns
        0 if a is -inf and b is not -inf;
        +inf if b is -inf (and a not -inf);
        exp (a) / exp (b) otherwise
    """
    if not a > float ("-inf") and not b > float ("-inf"):
        a = 0
        b = 0
    elif not b > float ("-inf"):
        return float ("+inf")
    
    answ = safe_exp (a - b)
    import math
    import traceback
    if math.isnan (answ):
        print ("LOLOLOLOLOLOLOLOL nan")
        print (a)
        print (b)
        print (traceback.print_stack ())
        print ("--------------")
    return answ


def safe_exp (a):
    """ Calculates np.exp (a) safely. 
    
    Parameters
        a: a float.

    Returns
        exp (a) if a is less than the log of the greatest float 
            possible in python, and +inf otherwise.
    """
    # python max float is
    max_float = 1.7976931348623177e+308
    # i'm not sure if that changes from machine to machine
    max_float = sys.float_info.max
    log_max = np.log (max_float)
    if a > log_max:
        return float ('inf')
    else:
        return np.exp (a)
    # ps: there's no warning when using super small arguments of np.exp


def clean_tag (xmlnode):
    """ Removes the namespace of an element tag. 
    
    Parameters
        xmlnode: an xml node object.

    Returns
        the node name, withouth the namespace tag.
    """
    return etree.QName (xmlnode.tag).localname


def create_dir_safe (dir_name):
    """ Safely creates a directory. 
    
    Parameters
        dir_name: a string with the name of the directory to be created.
    """
    try:
        os.makedirs(dir_name)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def plot_theta_var_sample (sample, var_idx, fig_name, plot_label, 
        custom_dir=None):
    """ Creates a plot of a sample of a multivariate variable.

    This function can plot an estimated distribution of a component of
    a multivariate random variable. The estimation is performed using 
    kernel density estimation. It is expected that the sampled variable
    is multivariate and that the user desires to estimate only one of
    its components per call to this function.
    
    Parameters
        sample: a list of points of a multivariate variable. Each
            point is also represented by a list.
        var_idx: the index of the component of interest in the random 
            variable. 
        fig_name: the filename of the output figure.
        plot_label: the title of the plot.
        custom_dir: the directory in which the output will be saved.
    """
    if custom_dir is not None:
        create_dir_safe (custom_dir)
        print ("Custom dir: " + custom_dir)
        fig_name = custom_dir + "/" + fig_name
        
    x = np.array ([theta[var_idx].value for theta in sample])
    sns_plot = sns.kdeplot (x, label=plot_label)
    fig = sns_plot.get_figure ()
    plt.legend ()
    fig.savefig (fig_name)
    plt.clf ()
