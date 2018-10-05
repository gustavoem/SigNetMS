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
        root = etree.Element ("dataset")
        root.set ("noise", "normal")
        root.set ("name", "Automatic generated dataset.")


        tree = etree.ElementTree (root)

        tree.write(file_name, encoding='utf-8', standalone=True, xml_declaration=True) 
