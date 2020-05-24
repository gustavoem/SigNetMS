class RandomParameterList:
    """ This class stores a list of random parameters. This set of
        random parameters can be considered as candidates of model
        parameters that were used for a biological experiment.
        Condiering that, one of the attributes of objects of this class,
        __experimental_error, represents the standard deviation of
        experimental error when assuming the correct model parameters
        are the ones present on this object.

        Attributes
            __param_list (list): a list of RandomParameter objects
            __experimental_error (RandomParameter): a random parameter
                that represents the standart deviation of experimental
                error associated to this set of parameters.
            __iterator: a simple iterator to iterate in model
                parameters.
    """

    def __init__ (self):
        """ Default constructor. """
        self.__param_list = []
        self.__experimental_error = None
        self.__iterator = None


    def __getitem__ (self, key):
        """ Defines indexed access. """
        return self.__param_list[key]
    

    def __iter__ (self):
        """ Iterator start. """
        self.__iterator = iter (self.__param_list)
        return self.__iterator


    def __next__ (self):
        """ Iterator step. """
        return next (self.__iterator)

    
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
        """ Sets the experimental error.

            Parameters:
                p: a RandomParameter object that represents the
                    experimental error random variable.
        """
        self.__param_list.append (p)
        self.__experimental_error = p
    

    def append (self, c):
        """ Adds a parameter to the end of the list.

            Parameters
                c: a RandomParameter to be added to this list
        """
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
        """ Returns the values of the parameters on the list.
            Note: does not include experimental error.
        """
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

    
    def get_p (self):
        """ Given the parameters and its distributions (priors), returns
            the value of the pdf of the joint distribution of parameters
            on the point represented by this object. """
        prob = 1
        for p in self.__param_list:
            prob *= p.get_p ()
        return prob


    def get_log_p (self):
        """ Returns the log of the prior joint probability of 
            this parameter list. """
        log_prob = 0
        for p in self.__param_list:
            log_prob += p.get_log_p ()
        return log_prob


    def get_size (self):
        """ Returns the number of elements. """
        return len (self.__param_list)
