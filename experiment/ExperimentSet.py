import sys
sys.path.insert (0, '..')

from lxml import etree
from experiment.Experiment import Experiment
from utils import clean_tag
import numpy as np

class ExperimentSet:
    """ This class stores a set of experiments. """

    def __init__ (self, filename=""):
        """ Default constructor. """
        self.__experiment_set = []
        if filename != "":
            self.load_data_file (filename)
        self.__iterator = None
    
    
    def __getitem__ (self, key):
        """ Defines indexed access. """
        return self.__experiment_set[key]


    def __iter__ (self):
        """ Iterator start. """
        self.__iterator = iter (self.__experiment_set)
        return self.__iterator


    def __next__ (self):
        """ Iterator step. """
        return next (self.__iterator)


    def add (self, experiment):
        """ Adds one experiment. """
        self.__experiment_set.append (experiment)


    def get_size (self):
        """ Returns the number of experiments in this set. """
        return len (self.__experiment_set)
    
    
    def load_data_file (self, file_name):
        """ This method reads an experiment set file and add all 
            experiments found to this object. """
        tree = etree.parse (file_name)
        root = tree.getroot ()
    
        if clean_tag (root) != "ExperimentSet":
            print ("Wrong experiment data syntax. Root tag should be" \
                + "<ExperimentSet>")
            return 

        experiments_arr = []
        for experiment_tag in root.getchildren ():
            rows = []

            if clean_tag (experiment_tag) != "Experiment":
                print ("Wrong experiment data syntax. The children of" \
                        + " <ExperimentSet> can only be of tag" \
                        + " <Experiment>.")
        
            for children in experiment_tag.getchildren ():
                if clean_tag (children) == "row":
                    row = self.__read_xml_row (children, file_name)
                    rows.append (row)
                elif clean_tag (children) == "condition" :
                    continue
                elif clean_tag (children) == "interpretation":
                    interp = self.__read_interpretation (children)
                else:
                    print ("Unexpected child of dataset in" + file_name)
            rows = np.array (rows)
    
            time_idx = interp.index ("time")
            times = rows[:, time_idx]
            for i in range (len (interp)):
                if i == time_idx:
                    continue
                expression = interp[i]
                var_values = rows[:, i]
                experiment = Experiment (times, var_values, expression)
                experiments_arr.append (experiment)

        for e in experiments_arr:
            self.add (e)

    @staticmethod
    def __read_interpretation (interp):
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
    

    @staticmethod
    def __read_xml_row (row_tag, file_name):
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

    
    def save_to_file (self, file_name):
        """ Saves experiment set into an XML file. """
        experiments = self.__experiment_set
        root = etree.Element ("ExperimentSet")
        
        for exp in experiments:
            exp_root = etree.SubElement (root, "Experiment")
            for i in range (len (exp.times)):
                row = etree.SubElement (exp_root, "row")
                time_elm = etree.SubElement (row, "element")
                time_elm.set ("value", str (exp.times[i]))
                time_elm.set ("index", "0")
                value_elm = etree.SubElement (row, "element")
                value_elm.set ("value", str (exp.values[i]))
                value_elm.set ("index", "1")
            etree.SubElement (exp_root, "condition")
            interp_elm = etree.SubElement (exp_root, "interpretation")
            time_interp = etree.SubElement (interp_elm, "time")
            time_interp.set ("col", "0")
            var_interp = etree.SubElement (interp_elm, "readout")
            var_interp.set ("expression", exp.measure_expression)
            var_interp.set ("col", "1")
            
        tree = etree.ElementTree (root)
        tree.write (file_name, pretty_print=True, encoding='utf-8',
                standalone=True, xml_declaration=True)

    def get_as_abcsysbio_syntax (self):
        """ Print each of the experiments to standard output using the
            ABC-SysBio syntax. """
        out = ''
        i = 1
        for exp in self.__experiment_set:
            out += exp.get_as_abcsysbio_syntax (i=i) + '\n'
            i += 1
        return out
