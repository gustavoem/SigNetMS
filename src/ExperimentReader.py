from Experiment import Experiment
from lxml import etree
from utils import clean_tag

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
    
    if clean_tag (root) != "dataset":
        print ("Wrong experiment data syntax. Root tag should be \
                <dataset>")
        return None
    
    for children in root.getchildren ():
        if clean_tag (children) == "data":
            values = read_xml_rows (children, file_name)
        elif clean_tag (children) == "condition" :
            continue
        elif clean_tag (children) == "interpretation":
            interpretation = read_interpretation (children, file_name)
        else:
            print ("Unexpected child of dataset in" + file_name)
    
    experiments = []
    time_idx = interpretation.index ("time")
    times = values[time_idx]
    for i in range (len (interpretation)):
        if i == time_idx:
            continue
        expression = interpretation[i]
        var_values = values[i]
        experiment = Experiment (times, var_values, expression)
        experiments.append (experiment)
    return experiments


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


def read_xml_rows (rows, file_name):
    """ Reads experiment rows on a data file. """
    columns = len (rows[0].getchildren ())
    values = [[] for x in range (columns)]

    for row in rows:
        for element in row.getchildren ():
            if clean_tag (element) != "element":
                print ("Unexpected child of row in " + file_name)
            attribs = element.attrib
            index = int (attribs["index"])
            value = float (attribs["value"])
            values[index].append (value)
    return values
