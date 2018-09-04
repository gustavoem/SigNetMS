class Experiment:
    """ This class stores experiments. """

    def __init__ (self, times, values, var):
        """ Default constructor. The list times represents the times
            of each read in values. var is the name of the variable read
            on those reads. """
        self.var = var
        self.times = times
        self.values = values
