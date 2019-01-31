import re
import matplotlib
matplotlib.use('Agg')

from scipy.integrate import odeint
from scipy.interpolate import spline
from sympy import diff
import matplotlib.pyplot as plt
import numpy as np


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


    def evaluate_on (self, time_points, initial_state_map=None):
        """ Returns the state of the systems variables at the specified
            time points. initial_state_map is an optional parameter 
            that shoud contain a map variable -> value that is used
            to modify the initial state of a few (or even all) variables
            of the system. """
        time_points = np.array (time_points)
        zeroed_times = False
        if time_points[0] != 0:
            time_points = np.insert (time_points, 0, 0)
            zeroed_times = True

        initial_state = self.initial_state.copy ()
        if initial_state_map != None:
            for var in initial_state_map:
                idx = self.index_map[var]
                initial_state[idx] = initial_state_map[var]

        sys_function = self.__create_system_function ()
        y, infodict = odeint (sys_function, initial_state, time_points, 
                mxstep=1000,full_output=True)

        values_map = {}
        for var in self.index_map:
            idx = self.index_map[var]
            if zeroed_times:
                # ignore first entry (initial state)
                values_map[var] = list (y[1:, idx])
            else:
                values_map[var] = list (y[:, idx])
        return values_map


    def overtime_plot (self, var_list, t, initial_state_map=None, 
            filename=''):
        """ Plots the values of a variable VAR over time t. """
        values_map = self.evaluate_on (t, initial_state_map)
        ODES.__plot (values_map, var_list, t, filename)


    # possible speedup: call this function only when it's necessary
    # to redefine the system
    def __create_system_function (self):
        """ Creates a function that describes the dynamics of the 
            system. """

        evaluable_formulas = []
        for formula in self.rate_eq:
            formula = ODES.__make_evaluable (formula)
            evaluable_formulas.append (formula)

        # Parameters are constant over time
        symbol_table = dict (self.param_table)

        def system_function (state, t):
            dstatedt = []
            current_state = {}

            for var, idx in self.index_map.items ():
                symbol_table[var] = state[idx]

            for i in range (len (state)):
                x = ODES.__calc_evaluable_func (evaluable_formulas[i], 
                        symbol_table)
                dstatedt.append (x)
            return dstatedt

        return system_function


    def get_system_jacobian (self):
        """ Creates the jacobian of the function that describes the 
            dynamics of the system. """

        #system_function = f (state) = (f_1 (state), ..., f_n (state))
        J_functions = []
        
        # For each component of f
        for i in range (len (self.rate_eq)):
            J_functions.append ([])
            f_i = self.rate_eq[i]
            
            J_functions[i] = [None for x in range (len (self.rate_eq))]

            # For each var we calculate  df_i/dx_j
            for var in self.index_map:
                j = self.index_map[var]
                dfdvar = ODES.__derivate (f_i, var)
                J_ij = ODES.__formula_to_lambda (dfdvar)
                J_functions[i][j] = J_ij
        
        # J is a matrix of lambdas
        # define a function here that returns the J matrix evaluated
        symbol_table = dict (self.param_table)
        def jacobian_function (state, t):
            for var, idx in self.index_map.items ():
                symbol_table[var] = state[idx]

            J_values = []
            for i in range (len (state)):
                J_values.append ([])
                for j in range (len (state)):
                    J_values[i].append (J_functions[i][j] (symbol_table))
            return J_values

        return jacobian_function


    @staticmethod 
    def __derivate (f, x):
        """ Returns the derivative of f in respect to x. This function
            only works for f that depends only linearly (or inverse 
            linearly) on x. """
        var_pattern = r'(^|\s|\+|-|\*|\/)(' + x + ')($|\s|\+|-|\*|\/)'
        var = x
        fX = re.sub (var_pattern, lambda m: m.group (1) + ' X ' + \
                m.group (3), f)
        dfX = str (diff (fX, 'X'))
        var_pattern = r'(^|\s|\+|-|\*|\/)(X)($|\s|\+|-|\*|\/)'
        df = re.sub (var_pattern, lambda m: m.group (1) + var + \
                m.group (3), dfX)
        return df


    @staticmethod
    def __formula_to_lambda (formula):
        """ Given a string formula, transform it on a lambda function. 
        """
        evaluable_formula = ODES.__make_evaluable (formula)
        def func (symbol_table):
            return ODES.__calc_evaluable_func (evaluable_formula, 
                    symbol_table)
        return func


    @staticmethod
    def __make_evaluable (formula):
        """ Transforms a formula in an evaluable formula. We do that by
            replacing a variable var by symbol_table[var] in the 
            formula. """
        new_formula = re.sub (r'(([A-z]|_)\w*)', 
                lambda m: 
                    m.group (0) if m.group (0) == "pow" 
                    else "symbol_table['" + m.group (0) + "']", formula)
        return new_formula


    @staticmethod
    def __calc_evaluable_func (func, symbol_table):
        """ Evaluate func in the scope of symbol table. """
        x = 0
        try:
            x = eval (func)
        except Exception as e:
            print ("Couldn't evaluate formula: \n " +
                    "\t" + func + "\n"
                    "Did you define system variables and " +
                    "parameters correctly?\n")
        return x


    @staticmethod
    def __plot (values_map, var_list, t, filename):
        """ Plots values of vars in var_list that were observed on time 
            t. """
        legend = []
        for var in var_list:
            legend.append (var)
            values = values_map[var]
            plt.plot (t, values)
        plt.legend (legend)
        if filename:
            plt.savefig (filename)
        else:
            plt.show ()
        plt.clf ()