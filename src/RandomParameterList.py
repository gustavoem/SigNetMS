import numpy as np

class RandomParameterList:
    """ This class stores a list of random parameters. """

    def __init__ (self):
        """ Default constructor. """
        self.__param_list = []
        self.__experimental_error = None


    def __getitem__ (self, key):
        """ Defines indexed access. """
        return self.__param_list[key]
    

    def __iter__ (self):
        """ Iterator start. """
        self.__current = 0
        return self


    def __next__ (self):
        """ Iterator step. """
        if self.__current >= len (self.__param_list):
            self.__current = 0
            raise StopIteration
        else:
            a = self.__param_list[self.__current]
            self.__current += 1
            return a

    
    def get_copy (self):
        """ Returns a copy of this parameter list. """
        copy = RandomParameterList ()
        last_idx = len (self.__param_list)

        if self.__experimental_error != None:
            last_idx = -1
            sigma_cpy = self.__experimental_error.copy ()
            copy.set_experimental_error (sigma_cpy)

        for p in self.__param_list[:last_idx]:
            copy.append (p.copy ()) 
        return copy


    def set_experimental_error (self, p):
        """ Sets the random parameter that represents the experimental
            error. """
        self.__param_list.append (p)
        self.__experimental_error = p
    

    def append (self, c):
        """ Adds a parameter to the end of the list. """
        if self.__experimental_error == None:
            insert_index = len (self.__param_list)
        else:
            insert_index = -1

        if type (c) == list:
            for p in c:
                self.__param_list.insert (insert_index, p)
        else:
            self.__param_list.insert (insert_index, c)


    def get_values (self):
        """ Returns the values of the parameters on the list. """
        values = []
        for p in self.__param_list:
            values.append (p.value)
        return values


    def get_experimental_error (self):
        """ Returns the experimental error value. """
        return self.__experimental_error.value


    def get_experimental_error_distribution (self):
        """ Returns a copy of the RandomParameter that defines the 
            distribution of the experimental error. """
        return self.__experimental_error.copy ()

    
    def get_model_parameters (self):
        """ Returns all but the experimental error parameters. """
        if self.__experimental_error == None:
            idx = len (self.__param_list)
        else:
            idx = len (self.__param_list) - 1
        return self.__param_list[:idx]


    def get_size (self):
        """ Returns the number of elements. """
        return len (self.__param_list)
