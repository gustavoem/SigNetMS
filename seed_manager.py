import numpy

__program_seed__ = None

def set_seed (seed):
    """ Sets numpy.random seed.
    
        Parameters
            seed: a number
    """
    global __program_seed__
    if __program_seed__ == None:
        numpy.random.seed (seed)
        __program_seed__ = seed

def get_seed ():
    """ Returns the set seed. """
    global __program_seed__
    return __program_seed__
