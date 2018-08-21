import libsbml

class SBML:
    """ This class contains a representation for SBML objects """

    def __init__ (self):
        """ Default constructor """
        self.sbml_obj = None

    def load_file (self, file_name):
        """ Given an xml file construct an sbml object """
        reader = libsbml.SBMLReader ()
        sbmldoc = reader.readSBML (file_name)
        
        # if sbmldoc.getNumErrors() > 0:
            # if sbmldoc.getError(0).getErrorId() == libsbml.XMLFileUnreadable:
                # sbmldoc.printErrors()
            # elif sbmldoc.getError(0).getErrorId() == libsbml.XMLFileOperationError:
                # sbmldoc.printErrors()
            # else:
                # sbmldoc.printErrors()
            # sys.exit(1)


    def get_name (self):
        """ Returns model name """
        return "Some model"

    # @staticmethod
    # def __check_sbml_parsing_err (sbmldoc):


