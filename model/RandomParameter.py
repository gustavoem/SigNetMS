class RandomParameter:
    """ This class stores a random parameter. 

        Attributes
            name (string): the name of this random parameter.
            value (double): the value of this parameter.
            __distribution (object): an object with the distribution of
                this random parameter.
    """
    
    def __init__ (self, name, distribution):
        """ Default constructor.
        
            Parameters
                name: the random parameter name.
                distribution: a distribution object.
        """
        self.name = name
        self.value = None
        self.__distribution = distribution
        self.set_rand_value ()


    def copy (self):
        """ Returns a copy of the random parameter. """
        distribution_cpy = self.__distribution.copy ()
        new_p = RandomParameter (self.name, distribution_cpy)
        new_p.value = self.value
        return new_p


    def set_rand_value (self):
        """ Sets a random value for the parameter.
            
            This value is set randomly according to the distribution of
            the random parameter.
        """
        self.value = self.__distribution.rvs ()
        return self.value


    def get_distribution (self):
        """ Returns the distribution object of this random parameter. 
        """
        return self.__distribution

    
    def get_p (self):
        """ Returns the likelihood of the current parameter value.
            
            This method returns the value of the probability density
            function of the parameter on its curent value.
        """
        return self.__distribution.pdf (self.value)


    def get_log_p (self):
        """ Returns the log of the likelihood of the current parameter
            value. 
        """
        return self.__distribution.log_pdf (self.value)

