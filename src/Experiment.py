from lxml import etree

class Experiment:
    """ This class stores experiments. """

    def __init__ (self, times, values, var):
        """ Default constructor. The list times represents the times
            of each read in values. var is the name of the variable read
            on those reads. """
        self.var = var
        self.times = times
        self.values = values


    def save_to_file (self, file_name):
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
        var_interp.set ("expression", self.var)
        var_interp.set ("col", "1")

        tree = etree.ElementTree (root)
        tree.write(file_name, pretty_print=True, encoding='utf-8', standalone=True, xml_declaration=True) 
