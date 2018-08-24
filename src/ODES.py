class ODES:
    """ This class contains a representation for systems of ordinary
        differential equations. """

    def __init__ (self):
        """ Default constructor. """
        self.index_map = {}
        self.rate_eq = {}
        self.initial_state = []

    
    def add_equation (self, var, string):
        """ Adds an equation representing the change rate of a variable.
        """
        if var in self.index_map:
            self.rate_eq[var] = string
            return
        
        idx = len (self.index_map)
        self.index_map[var] = idx
        self.initial_state.append (None)


    def define_initial_value (self, var, value):
        """ Defines the initial value of a variable. """
        idx = self.index_map[var]
        self.initial_state[idx] = value


    def evaluate_on (self, time_points):
        """ Returns the state of the systems variables at the specified
            time points. """
        return self.initial_state


