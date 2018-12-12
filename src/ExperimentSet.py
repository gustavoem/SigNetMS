from lxml import etree
from Experiment import Experiment

class ExperimentSet:
    """ This class stores a set of experiments. """

    def __init__ (self):
        """ Default constructor. """
        self.__experiment_set = []
    
    
    def __iter__ (self):
        """ Iterator start. """
        self.__current = 0
        return self


    def __next__ (self):
        """ Iterator step. """
        if self.__current >= len (self.__experiment_set):
            self.__current = 0
            raise StopIteration
        else:
            a = sel.__experiment_set[self.__current]
            self.__current += 1
            return a


    def add (self, experiment):
        """ Adds one experiment. """
        self.__experiment_set.append (experiment)

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
            condition_elm = etree.SubElement (exp_root, "condition")
            interp_elm = etree.SubElement (exp_root, "interpretation")
            time_interp = etree.SubElement (interp_elm, "time")
            time_interp.set ("col", "0")
            var_interp = etree.SubElement (interp_elm, "readout")
            var_interp.set ("expression", exp.measure_expression)
            var_interp.set ("col", "1")
            
        tree = etree.ElementTree (root)
        tree.write (file_name, pretty_print=True, encoding='utf-8',
                standalone=True, xml_declaration=True)
