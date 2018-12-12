from Experiment import Experiment
from ExperimentSet import ExperimentSet
from lxml import etree
from utils import clean_tag
import numpy as np

def read_txt_experiment_file (file_name, var):
    """ This function can read experiment files. The experiment data is 
        organized in a list of Experiment objects. """
    f = open (file_name)
    experiment_values = []
    times = []
    experiments = []
    for line in f:
        l = line.split ()
        # remove comments
        if '#' in l:
            l = l[0:l.index ('#')]
        if len (l) == 0:
            # create experiment objects
            for exp_val in experiment_values:
                exp = Experiment (times, exp_val, var)
                experiments.append (exp)
                experiment_values = []
            times = []
            continue
        
        times.append (l[0])
        if (experiment_values == []):
            nof_exp_per_line = len (l) - 1
            experiment_values = [[] for i in range (len (l) - 1)]
        for i in range (1, len (l)):
            experiment_values[i - 1].append (float (l[i]))

    for exp_val in experiment_values:
        exp = Experiment (times, exp_val, var)
        experiments.append (exp)
        experiment_values = []

    f.close ()
    return experiments


def read_data_experiment_file (file_name):
    """ This function can read data experiment files. The experiment 
        data is organized in a list of Experiment objects. """
    tree = etree.parse (file_name)
    root = tree.getroot ()
    
    if clean_tag (root) != "ExperimentSet":
        print ("Wrong experiment data syntax. Root tag should be" \
                + "<ExperimentSet>")
        return None

    experiment_set = ExperimentSet () 
    for experiment_tag in root.getchildren ():
        rows = []

        if clean_tag (experiment_tag) != "Experiment":
            print ("Wrong experiment data syntax. The children of" \
                    + " <ExperimentSet> can only be of tag" \
                    + " <Experiment>.")
        
        for children in experiment_tag.getchildren ():
            if clean_tag (children) == "row":
                rows.append (read_xml_row (children, file_name))
            elif clean_tag (children) == "condition" :
                continue
            elif clean_tag (children) == "interpretation":
                interpretation = read_interpretation (children, 
                        file_name)
            else:
                print ("Unexpected child of dataset in" + file_name)
        rows = np.array (rows)
    
        time_idx = interpretation.index ("time")
        times = rows[:, time_idx]
        for i in range (len (interpretation)):
            if i == time_idx:
                continue
            expression = interpretation[i]
            var_values = rows[:, i]
            experiment = Experiment (times, var_values, expression)
            experiment_set.add (experiment)
    return experiment_set


def read_interpretation (interp, file_name):
    """ Reads the interpretation subtree of an experiment data file.  
    """ 
    interp_arr = [None] * len (interp)
    for element in interp:
        if (clean_tag (element) == "time"):
            idx = int (element.attrib["col"])
            interp_arr[idx] = "time"
        else:
            idx = int (element.attrib["col"])
            label = element.attrib["expression"]
            interp_arr[idx] = label
    return interp_arr


def read_xml_row (row_tag, file_name):
    """ Reads experiment rows on a data file. """
    columns = len (row_tag.getchildren ())
    row = [None for x in range (columns)]
    for element in row_tag.getchildren ():
        if clean_tag (element) != "element":
            print ("Unexpected child of row in " + file_name)
        attribs = element.attrib
        index = int (attribs["index"])
        value = float (attribs["value"])
        row[index] = value
    return row 
