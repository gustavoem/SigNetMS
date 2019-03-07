import matplotlib
matplotlib.use('Agg')

from lxml import etree
import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def safe_power (a, b):
    """ Calculates a ** b safely."""
    try:
        x = a ** b
    except OverflowError:
        x = float ("inf")
    return x


def safe_exp (a):
    """ Calculates np.exp (a) safely. """
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
    """ Removes the namespace of an element tag. """
    return etree.QName (xmlnode.tag).localname


def plot_theta_var_sample (sample, var_idx, fig_name, plot_label):
    x = np.array ([theta[var_idx].value for theta in sample])
    sns_plot = sns.kdeplot (x, label=plot_label)
    fig = sns_plot.get_figure ()
    plt.legend ()
    fig.savefig (fig_name)
    plt.clf ()
