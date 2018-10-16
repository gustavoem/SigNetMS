import numpy as np

class RandomParameterList:
    """ This class stores a list of random parameters. """

    def __init__ (self):
        """ Default constructor. """
        self.__param_list = []


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
        for p in self.__param_list:
            copy.append (p.copy ())
        return copy


    def append (self, c):
        """ Adds a parameter to the end of the list. """
        if type (c) == list:
            for p in c:
                self.__param_list.append (p)
        else:
            self.__param_list.append (c)


    def get_values (self):
        """ Returns the values of the parameters on the list. """
        values = []
        for p in self.__param_list:
            values.append (p.value)
        return values


    def get_size (self):
        """ Returns the number of elements. """
        return len (self.__param_list)
