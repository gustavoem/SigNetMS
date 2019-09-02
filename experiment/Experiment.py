from lxml import etree

class Experiment:
    """ This class represent an experiment. """

    def __init__ (self, times, values, measure):
        """ Default constructor. 
        
        Parameters
            times: a list that represents the time moments of each read 
                in the experiment.
            values: the value measured in each moment.
            measure: an expression that represents what was the 
                measurement for the experiment. This measure is 
                generally a mathematical expression written in terms of
                concentrations of chemical species.
        """
        self.measure_expression = measure
        self.times = times
        self.values = values
    

    def save_as_BioBayes_file (self, filename):
        """ Creates an xml file with the experiment data. This file
            is compatible with BioBayes software. """
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
        
        etree.SubElement (root, "condition")
        interpretation = etree.SubElement (root, "interpretation")
        time_interp = etree.SubElement (interpretation, "time")
        time_interp.set ("col", "0")
        var_interp = etree.SubElement (interpretation, "readout")
        var_interp.set ("expression", self.measure_expression)
        var_interp.set ("col", "1")

        tree = etree.ElementTree (root)
        tree.write(filename, pretty_print=True, encoding='utf-8', 
                standalone=True, xml_declaration=True)


    def get_as_abcsysbio_syntax (self, i=None):
        """ Prints the experiment to standard output using the 
            ABC-SysBio syntax. """
        if i is None:
            i = 1
        
        i_str = str (i)
        out = '<var' + i_str + '> '
        for val in self.values:
            out += str (val) + ' '
        out += '</var' + i_str + '>' 
        return out
