import numpy as np

class Constant:
    """ This class defines a constant distribution. """

    def __init__ (self, constant):
        self.__constant = constant


    def copy (self):
        cpy = Constant (self.__constant)
        return cpy


    def get_params (self):
        params = {}
        params["constant"] = self.__constant
        return params


    def mean (self):
        return self.__constant


    def rvs (self, n=None):
        if n is None:
            return self.__constant
        else:
            return np.array([self.__constant for _ in range (n)])


    def pdf (self, x):
        if abs (x - self.__constant) < 1e-4:
            return 1
        return 0


    def log_pdf (self, x):
        if self.pdf(x):
            return 0
        else:
            return float ("-inf")
