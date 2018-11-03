import sys
sys.path.insert (0, '../src/')

import unittest
from PriorsReader import read_priors_file

class TestPriorsReader (unittest.TestCase):

    def test_reader (self):
        """ Tests if the module can correctly read a prior definition 
            file. """
        priors = read_priors_file ('input/simple_enzymatic.priors')
        for x in priors:
            if x.name == 'k1':
                self.assertEqual (x.get_a (), 2.0)
                self.assertEqual (x.get_b (), 0.01)
            elif x.name == 'd1' or x.name == 'kcat':
                self.assertEqual (x.get_a (), 2.0)
                self.assertEqual (x.get_b (), 0.1)
            else:
                self.fail ()

    def test_sigma_prior (self):
        """ Tests if the module can read the experimental error 
            prior. """
        priors = read_priors_file ('input/simple_enzymatic.priors')
        sigma = priors.get_experimental_error ()
        self.assertEqual (sigma.get_a (), 2.0)
        self.assertEqual (sigma.get_b (), 2.6)
