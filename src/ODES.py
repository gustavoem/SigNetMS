import re

from scipy.integrate import odeint

class ODES:
    """ This class contains a representation for systems of ordinary
        differential equations. """

    def __init__ (self):
        """ Default constructor. """
        # A map var -> var index
        self.index_map = {}
        
        # The rate equation of a var
        self.rate_eq = []
        
        # The initial state of a var
        self.initial_state = []
    
        # A hash table with values of parameters
        self.param_table = {}


    def __create_var (self, var):
        """ Creates a new variable and adds it to the system. """
        idx = len (self.index_map)
        self.index_map[var] = idx

    
    def add_equation (self, var, formula):
        """ Adds an equation representing the change rate of a variable.
        """
        if var in self.index_map:
            self.rate_eq[index_map[var]] = formula
            return
        
        idx = self.__create_var (var)
        self.rate_eq.append (formula)
        self.initial_state.append (None)


    def define_initial_value (self, var, value):
        """ Defines the initial value of a variable. """
        idx = None
        if var not in self.index_map:
            idx = __create_var (var)
        else:
            idx = self.index_map[var]
        self.initial_state[idx] = value


    def define_parameter (self, param, value):
        """ Defines the value of some parameter. """
        self.param_table[param] = value


    def get_all_parameters (self):
        """ Returns all parameters of the system. """
        return self.param_table


    def evaluate_on (self, time_points):
        """ Returns the state of the systems variables at the specified
            time points. """
        sys_function = self.__create_system_function ()
        y = odeint (sys_function, self.initial_state, time_points)
        values_map = {}
        for var in self.index_map:
            idx = self.index_map[var]
            values_map[var] = list (y[:,idx])
        return values_map


    # possible speedup: call this function only when it's necessary
    # to redefine the system
    def __create_system_function (self):
        """ Creates a function that describes the dynamics of the 
            system. """

        evaluable_formulas = []
        for formula in self.rate_eq:
            formula = re.sub (r'(([A-z]|_)\w*)', 
                    lambda m: "symbol_table['" + m.group (0) + "']", 
                    formula)
            evaluable_formulas.append (formula)

        # Param state is constant over time
        symbol_table = dict (self.param_table)


        def system_function (state, t):
            dstatedt = []
            current_state = {}

            for var in self.index_map:
                idx = self.index_map[var]
                symbol_table[var] = state[idx]

            for i in range (len (state)):
                x = 0
                try:
                    x = eval (evaluable_formulas[i])
                except:
                    print ("Couldn't evaluate system rate formula. " +
                           "Did you define system variables and " +
                           "parameters correctly?")
                dstatedt.append (x)
            return dstatedt

        return system_function
