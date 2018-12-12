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


    def save_to_file (self, filename):
        """ Creates an xml file with the experiment data. """
        DATA_NAMESPACE = "http://bisb.gla.org/dataset"
        NSMAP = {None: DATA_NAMESPACE}

        root = etree.Element ("dataset", nsmap=NSMAP)
        root.set ("noise", "normal")
        root.set ("name", "Automatic generated dataset.")

        data = etree.SubElement (root, "data")
        data.set ("cols", "2")
        for i in range (len (self.times)):
            row = etree.SubElement (data, "row")
            element0 = etree.SubElement (row, "element")
            element0.set ("value", str (self.times[i]))
            element0.set ("index", "0")
            element1 = etree.SubElement (row, "element")
            element1.set ("value", str (self.values[i]))
            element1.set ("index", "1")
        
        condition = etree.SubElement (root, "condition")
        interpretation = etree.SubElement (root, "interpretation")
        time_interp = etree.SubElement (interpretation, "time")
        time_interp.set ("col", "0")
        var_interp = etree.SubElement (interpretation, "readout")
        var_interp.set ("expression", self.measure_expression)
        var_interp.set ("col", "1")

        tree = etree.ElementTree (root)
        tree.write(file_name, pretty_print=True, encoding='utf-8', 
                standalone=True, xml_declaration=True) 
