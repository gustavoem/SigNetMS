from lxml import etree

class Experiment:
    """ This class stores experiments. """

    def __init__ (self, times, values, measure):
        """ Default constructor. The list times represents the times
            of each read in values. var is the name of the variable read
            on those reads. """
        self.measure_expression = measure
        self.times = times
        self.values = values
