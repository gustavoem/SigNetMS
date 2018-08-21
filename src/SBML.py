import libsbml

class SBML:
    """ This class contains a representation for SBML objects. """

    def __init__ (self):
        """ Default constructor. """
        self.sbml_obj = None

    def load_file (self, file_name):
        """ Given an xml file construct an sbml object. """
        reader = libsbml.SBMLReader ()
        sbmldoc = reader.readSBML (file_name)
        
            
    def get_species_list (self):
        """ Returns a list with all species. """
        return []

    def get_name (self):
        """ Returns model name. """
        return "Some model"

    # @staticmethod
    # def __check_sbml_parsing_err (sbmldoc):

