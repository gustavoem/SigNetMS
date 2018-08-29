import re

from scipy.integrate import odeint

class ODES:
    """ This class contains a representation for systems of ordinary
        differential equations. """

    def __init__ (self):
        """ Default constructor. """
        self.index_map = {}
        self.rate_eq = []
        self.initial_state = []
        self.param_table = {}


    # def __create_var (self):

    
    def add_equation (self, var, formula):
        """ Adds an equation representing the change rate of a variable.
        """
        if var in self.index_map:
            self.rate_eq[index_map[var]] = formula
            return
        
        idx = len (self.index_map)
        self.index_map[var] = idx
        self.rate_eq.append (formula)
        self.initial_state.append (None)


    def define_initial_value (self, var, value):
        """ Defines the initial value of a variable. """
        # if var not in self.inex_map:
            # idx = len ()

        idx = self.index_map[var]
        self.initial_state[idx] = value


    def define_parameter (self, param, value):
        """ Defines the value of some parameter """
        pass


    def evaluate_on (self, time_points):
        """ Returns the state of the systems variables at the specified
            time points. """
        sys_function = self.__create_system_function ()
        y = odeint (sys_function, self.initial_state, time_points)
        return y

    
    # possible speedup: call this function only when it's necessary
    # to redefine the system
    def __create_system_function (self):
        """ Creates a function that describes the dynamics of the 
            system. """

        evaluable_formulas = []
        i = 0
        for formula in self.rate_eq:
            formula = re.sub (r'(([A-z]|_)\w*)', 
                    lambda m: "current_state['" + m.group (0) + "']", 
                    formula)
            evaluable_formulas.append (formula)
            i += 1


        def system_function (state, t):
            dstatedt = []
            current_state = {}

            for var in self.index_map:
                idx = self.index_map[var]
                current_state[var] = state[idx]

            for i in range (len (state)):
                x = 0
                try:
                    x = eval (evaluable_formulas[i])
                except:
                    print ("Couldn't evaluate system rate formula. " +
                           "Did you define system variables " + 
                           "correctly?")
                dstatedt.append (x)
            return dstatedt

        return system_function
