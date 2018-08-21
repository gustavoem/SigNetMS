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
        if self.__check_sbml_parsing_err (sbmldoc):
            return False
        self.sbml_obj = sbmldoc
        return True
        
            
    def get_species_list (self):
        """ Returns a list with all species. """
        model = self.sbml_obj.model
        number_of_species = model.getNumSpecies ()
        
        species = []
        for i in range (number_of_species):
            species.append (model.getSpecies (i).getName ())

        return species


    def get_name (self):
        """ Returns model name. """
        return self.sbml_obj.model.getName ()


    @staticmethod
    def __check_sbml_parsing_err (sbmldoc):
        """ Verifies if an SBMLReader () was succesful """
        if sbmldoc.getNumErrors () > 0:
            sbmldoc.printErrors ()
            return True
        return False
