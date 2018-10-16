import numpy as np

class RandomParameterList:
    """ This class stores a list of random parameters. """

    def __init__ (self):
        """ Default constructor. """
        pass
        self.__param_list = []

    
    def get_copy (self):
        """ Returns a copy of this parameter list. """
        # copy = RandomParameterList ()
        # for p in self.__param_list
        pass


    def append (self, p):
        self.__param_list.append (p)


    def get_params_values (self):
        return None


    def get_size (self):
        """ Returns the number of elements. """
        return len (self.__param_list)
