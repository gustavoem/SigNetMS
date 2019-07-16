import re
import matplotlib
matplotlib.use('Agg') #uncomment when running on server

from scipy.integrate import odeint
import scipy.integrate as spi
import sympy as sym
from sympy.parsing.sympy_parser import parse_expr
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

        # A sympy object that represents the system
        self.sys_eq = None

        # A sympy object that represents the system variables in sys_eq
        self.sys_vars = None

        # A sympy object that represents the system parameters in sys_eq
        self.sys_params = None


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
        self.sys_eq = None
        self.sys_vars = None
        self.sys_params = None


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
        self.sys_eq = None
        self.sys_vars = None
        self.sys_params = None


    def get_all_parameters (self):
        """ Returns all parameters of the system. """
        return self.param_table


    def __integrate_with_odeint (self, sys_f, initial_state, 
            time_points):
        """ Integrates using scipy odeint """
        # The order of iterations of a dict in python will not change
        # if nothing has been added (if that is true, the system is 
        # rewritten).
        args = [self.param_table[param] for param in self.param_table]
        y, _ = odeint (sys_f, initial_state, time_points, args=(args,),
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


    def __get_sym_vars_equations (self):
        """ Returns a tuple with sympy symbols of vars, parameters and 
            its correspondant reaction rate equations. """
        local_dict = {}
        var_symbols = []
        for var in self.index_map:
            var_sym = sym.symbols (var)
            var_symbols.append (var_sym)
            local_dict[var] = var_sym

        parameters = []
        for param in self.param_table:
            p_symbol = sym.symbols (param)
            parameters.append (p_symbol)
            local_dict[param] = p_symbol

        rate_equations = []
        for var in self.index_map:
            var_idx = self.index_map[var]
            expr = parse_expr (self.rate_eq[var_idx], local_dict=local_dict)
            rate_equations.append (expr)
        return var_symbols, parameters, rate_equations


    def __define_sys_eq (self):
        """ Creates a sympy object that represents the system of 
            ordinary differential equations. """
        n = len (self.rate_eq)
        m = len (self.param_table)
        var_names, params, equations = self.__get_sym_vars_equations ()
        brute_odes_rhs = sym.Matrix (equations) # var, param not indexed
        sym_vars = sym.MatrixSymbol ('y', n, 1)
        sym_dvars = sym.MatrixSymbol ('dy', n, 1)
        sym_p = sym.MatrixSymbol ('p', m, 1)
        var_to_sym_map = dict (zip (var_names, sym_vars))
        p_to_sym_map = dict (zip (params, sym_p))
        odes_rhs_idxdvar = brute_odes_rhs.xreplace (var_to_sym_map)
        odes_rhs = odes_rhs_idxdvar.xreplace (p_to_sym_map)
        fun = sym.lambdify ([sym_vars, sym_p], odes_rhs)
        self.sys_eq = sym.Eq (sym_dvars, odes_rhs)
        self.sys_vars = sym_vars
        self.sys_params = sym_p


    def __create_system_function (self):
        """ Creates a function that describes the dynamics of the 
            system. """
        if self.sys_eq == None:
            self.__define_sys_eq ()
        
        sys_eq_rhs = self.sys_eq.rhs
        sys_vars = self.sys_vars
        sys_params = self.sys_params
        sys_lambda = sym.lambdify ([sys_vars, sys_params], sys_eq_rhs)
        # sys_lambda is a lambda function that was automatically 
        # generated by sympy. It takes as arguments a numpy array of
        # shape (n, 1) (states) and another with shape (m, 1) 
        # (parameters). It returns also a numpy array, corresponding to
        # the evaluated derivatives, of shape (n, 1).
        wrapped_sys_fun = self.odeint_sys_wrapper (sys_lambda)
        return wrapped_sys_fun 


    def odeint_sys_wrapper (self, lamb):
        """ This is a wrapper to the lambda functions that were created
            with sympy to represent the system function and also its 
            Jacobian. 
            
            The lambda function is expected to receive as input a numpy 
            array of shape (n, 1) (states) and another with shape (m, 1) 
            (parameters); it should also return a numpy array, 
            corresponding to the evaluated derivatives (or second 
            derivatives in the case of the Jacobian), of shape (n, 1). 
            
            The wrapper function is expected to receive two arguments
            """
        n = len (self.rate_eq)
        m = len (self.param_table)
        if m == 0:
            def wrapped_fun (t, state, args):
                npstate = np.array (state)
                npstate.shape = (n, 1)
                answ = lamb (npstate, None)
                return answ.squeeze ()
        else:
            def wrapped_fun (t, state, args):
                npstate = np.array (state)
                npparams = np.array (args)
                npstate.shape = (n, 1)
                npparams.shape = (m, 1)
                answ = lamb (npstate, npparams)
                return answ.squeeze ()
        return wrapped_fun


    def get_system_jacobian (self):
        """ Creates the jacobian of the function that describes the 
            dynamics of the system. """
        # system_function = f (state) = (f_1 (state), ..., f_n (state))
        if self.sys_eq == None:
            self.__create_system_function () 

        rhs = self.sys_eq.rhs
        sys_vars = self.sys_vars
        sys_params = self.sys_params
        sym_jacobian = rhs.jacobian (sys_vars)
        lam_jacobian = sym.lambdify ([sys_vars, sys_params], \
                sym_jacobian)
        jac_fun = self.odeint_sys_wrapper (lam_jacobian)
        return jac_fun


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
