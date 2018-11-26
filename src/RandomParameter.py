class RandomParameter:
    """ This class stores a random parameter. """
    
    def __init__ (self, name, distribution):
        """ Default constructor. Instantiate a gamma distributed 
            parameter with shape a and scale parameter b. """
        self.name = name
        self.value = None
        self.__distribution = distribution
        self.set_rand_value ()


    def copy (self):
        """ Returns a copy of the self. """
        distribution_cpy = distribution.copy ()
        new_p = RandomParameter (self.name, distribution_cpy)
        new_p.value = self.value
        return new_p


    def set_rand_value (self):
        """ Sets the parameter a random value distributed as a 
            Gamma (a, b)."""
        self.value = distribution.rvs ()
        return self.value

    
    def get_p (self):
        """ Returns the probability of the current value, given the 
            distribution. """
        return self.__distribution.pdf (self.value)


    def get_a (self):
        """ Returns the shape of the gamma distribution of this random
            parameter. """
        return self.__a
    
    
    def get_b (self):
        """ Returns the inverse scale (shape) of the gamma distribution 
            of this random parameter. """
        return self.__b
