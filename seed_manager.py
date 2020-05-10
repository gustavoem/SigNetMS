import numpy

@singleton # pylint: disable=undefined-variable
class SeedManager:
    def __init__ (self):
        self.__seed = None

    def set_seed (self, seed):
        if self.__seed == None:
            numpy.random.seed (seed)
            self.__seed = seed

    def get_seed (self):
        return self.__seed


def set_seed (seed):
    """ Sets numpy.random seed.
    
        Parameters
            seed: a number
    """
    sm = SeedManager.instance ()
    sm.set_seed (seed)

def get_seed ():
    """ Returns the set seed. """
    return SeedManager.instance ().get_seed ()
