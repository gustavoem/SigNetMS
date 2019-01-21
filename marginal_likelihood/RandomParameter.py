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
        distribution_cpy = self.__distribution.copy ()
        new_p = RandomParameter (self.name, distribution_cpy)
        new_p.value = self.value
        return new_p


    def set_rand_value (self):
        """ Sets the parameter a random value distributed as a 
            Gamma (a, b)."""
        self.value = self.__distribution.rvs ()
        return self.value


    def get_distribution (self):
        """ Returns the distribution of this object. """
        return self.__distribution

    
    def get_p (self):
        """ Returns the probability of the current value, given the 
            distribution. """
        return self.__distribution.pdf (self.value)
