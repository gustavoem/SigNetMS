import re
import matplotlib
matplotlib.use('Agg') #uncomment when running on server

from scipy.integrate import odeint
import scipy.integrate as spi
from sympy import diff
from asteval import Interpreter
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

        # The model's name
        self.name = ''


    def __create_var (self, var):
        """ Creates a new variable and adds it to the system. """
        idx = len (self.index_map)
        self.index_map[var] = idx


    def print_equations (self):
        """ Prints all var equations. """
        for var in self.index_map:
            print ("d" + var + "/dt = ", end="")
            print (self.rate_eq[self.index_map[var]], end='\n\n')


    def add_equation (self, var, formula):
        """ Adds an equation representing the change rate of a variable.
        """
        if var in self.index_map:
            self.rate_eq[self.index_map[var]] = formula
            return
        
        self.__create_var (var)
        self.rate_eq.append (formula)
        self.initial_state.append (None)


    def define_initial_value (self, var, value):
        """ Defines the initial value of a variable. """
        idx = None
        if var not in self.index_map:
            idx = self.__create_var (var)
        else:
            idx = self.index_map[var]
        self.initial_state[idx] = value


    def define_parameter (self, param, value):
        """ Defines the value of some parameter. """
        self.param_table[param] = value


    def get_all_parameters (self):
        """ Returns all parameters of the system. """
        return self.param_table


    def __integrate_with_odeint (self, sys_f, initial_state, 
            time_points):
        """ Integrates using scipy odeint """
        y, _ = odeint (sys_f, initial_state, time_points, 
                mxstep=5000, full_output=True, tfirst=True, atol=1e-6,
                rtol=1e-8)
        return y


    def __integrate_with_stiff_alg (self, sys_f, initial_state, 
            time_points):
        """ Integrates using scipy.ode.integrate with an algorithm for
            stiff problems. """
        dt = .1
        ode = spi.ode (sys_f)
        ode.set_integrator ('vode', nsteps=5000, method='bdf', 
                max_step=dt, atol=1e-6, rtol=1e-4)
        ode.set_initial_value (initial_state, time_points[0])
        y = [initial_state]
        i = 1
        while ode.t < time_points[-1]: # and ode.succesful ():
            next_point = min (ode.t + dt, time_points[i])
            ode.integrate (next_point)
            if ode.t >= time_points[i]:
                y.append (ode.y)
                i += 1
        return np.array (y)


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
        y = self.__integrate_with_odeint (sys_function, 
                initial_state, time_points)
        
        values_map = {}
        for var in self.index_map:
            idx = self.index_map[var]
            if zeroed_times:
                # ignore first entry (initial state)
                values_map[var] = list (y[1:, idx])
            else:
                values_map[var] = list (y[:, idx])
        return values_map


    def evaluate_exp_on (self, exp, time_points, 
            initial_state_map=None):
        """ Integrates the system and returns an array containing the
            values of exp on each time-step evaluated of the system."""
        system_states = self.evaluate_on (time_points, \
                initial_state_map)
        aeval = Interpreter ()
        values = []
        for i in range (len (time_points)):
            for var in system_states:
                aeval.symtable[var] = system_states[var][i]
            value = aeval (exp)
            values.append (value)
        return values



    def overtime_plot (self, var_list, t, initial_state_map=None, 
            filename='', xlabel=None, ylabel=None, title=None):
        """ Plots the values of a variable VAR over time t. """
        values_map = self.evaluate_on (t, initial_state_map)
        ODES.__plot (values_map, var_list, t, filename, xlabel, ylabel,
                title)


    # possible speedup: call this function only when it's necessary
    # to redefine the system
    def __create_system_function (self):
        """ Creates a function that describes the dynamics of the 
            system. """

        # Start the symbol table of intepreter with parameter values
        # because they are constant over time.
        rate_eq_evaluator = Interpreter ()
        for param in self.param_table:
            rate_eq_evaluator.symtable[param] = self.param_table[param]
        
        #pylint: disable=unused-argument
        def system_function (t, state):
            dstatedt = []

            # Add variables states to the interpreter symbol table
            for var, idx in self.index_map.items ():
                rate_eq_evaluator.symtable[var] = state[idx]

            for i in range (len (state)):
                formula = self.rate_eq[i]
                x = ODES.__calc_func (formula, {}, rate_eq_evaluator)
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

        #pylint: disable=unused-argument
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
        var_pattern = r"(^|\s|\+|-|\*|\/)(" + x + r")($|\s|\+|-|\*|\/)"
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
        def func (symbol_table):
            return ODES.__calc_func (formula, symbol_table)
        return func


    @staticmethod
    def __calc_func (func, symbol_table, interpreter=None):
        """ Evaluate func in the scope of symbol table. """
        x = 0
        if interpreter == None:
            interpreter = Interpreter ()

        for symbol in symbol_table:
            interpreter.symtable[symbol] = symbol_table[symbol]
        
        try:
            x = interpreter (func)
        except NameError as e:
            print ("Couldn't evaluate formula: \n " +
                    "\t" + func + "\n"
                    "Did you define system variables and " +
                    "parameters correctly?\n")
            print (e)
        except SyntaxError as e:
            print ("Couldn't evaluate formula: \n " +
                    "\t" + func + "\n")
            print (e)
        return x


    @staticmethod
    def __plot (values_map, exp_list, t, filename, xlabel, ylabel, 
            title):
        """ Plots values of vars in var_list that were observed on time 
            t. """
        legend = []
        var_list = [var for var in values_map.keys ()]

        for exp in exp_list:
            legend.append (exp)
            values = []
            for i in range (len (t)):
                symbol_table = {}
                for var in var_list:
                    symbol_table[var] = values_map[var][i]
                v = ODES.__calc_func (exp, symbol_table)
                values.append (v)
            plt.plot (t, values)
        plt.legend (legend)
        
        if xlabel: 
            plt.xlabel (xlabel)
        if ylabel: 
            plt.ylabel (ylabel)
        if title: 
            plt.title (title)

        if filename:
            plt.savefig (filename)
        else:
            plt.show ()
        plt.clf ()
