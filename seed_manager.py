import numpy

# pylint: disable=undefined-variable
__program_seed__ = None

def set_seed (seed):
    """ Sets numpy.random seed.
    
        Parameters
            seed: a number
    """
    # pylint: disable=global-statement
    global __program_seed__
    if __program_seed__ == None:
        numpy.random.seed (seed)
        __program_seed__ = seed

def get_seed ():
    """ Returns the set seed. """
    # pylint: disable=global-statement
    global __program_seed__
    return __program_seed__
