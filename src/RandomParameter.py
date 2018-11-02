import numpy as np

class RandomParameter:
    """ This class stores a random parameter. """
    
    def __init__ (self, name, a, b):
        """ Default constructor. Instantiate a gamma distributed 
            parameter with shape a and scale parameter b. """
        self.name = name
        self.__a = a
        self.__b = b
        self.value = None
        self.set_rand_value ()


    def copy (self):
        """ Returns a copy of the self. """
        new_p = RandomParameter (self.name, self.__a, self.__b)
        new_p.value = self.value
        return new_p


    def set_rand_value (self):
        """ Sets the parameter a random value distributed as a 
            Gamma (a, b)."""
        self.value = np.random.gamma(self.__a, self.__b)
        return self.value


    def get_a (self):
        """ Returns the shape of the gamma distribution of this random
            parameter. """
        return self.__a
    
    
    def get_b (self):
        """ Returns the inverse scale (shape) of the gamma distribution 
            of this random parameter. """
        return self.__b
