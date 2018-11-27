import matplotlib
matplotlib.use('Agg')

from lxml import etree
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


def clean_tag (xmlnode):
    """ Removes the namespace of an element tag. """
    return etree.QName (xmlnode.tag).localname


def plot_theta_var_sample (sample, var_idx, fig_name):
    x = np.array ([theta[var_idx].value for theta in sample])
    sns_plot = sns.kdeplot (x)
    fig = sns_plot.get_figure ()
    fig.savefig (fig_name)
