import numpy as np

class RandomParameter:
    """ This class stores a random parameter. """
    
    def __init__ (self, name, a, b):
        """ Default constructor. Instantiate a gamma distributed 
            parameter with shape a and inverse scale parameter b. """
        self.name = name
        self.__a = a
        self.__b = b
        self.value = None
        self.set_rand_value ()


    def set_rand_value (self):
        """ Sets the parameter a random value distributed as a 
            Gamma (a, b)."""
        self.value = np.random.gamma(self.__a, 1.0 / self.__b)
        return self.value


    def get_a (self):
        """ Returns the shape of the gamma distribution of this random
            parameter. """
        return self.__a
    
    
    def get_b (self):
        """ Returns the scale of the gamma distribution of this random
            parameter. """
        return self.__b
