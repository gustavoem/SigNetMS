import sys
sys.path.insert (0, '..')

import unittest
from marginal_likelihood.PriorsReader import read_priors_file
from marginal_likelihood.PriorsReader import define_sbml_params_priors 
from marginal_likelihood.SBML import SBML

class TestPriorsReader (unittest.TestCase):

    def test_reader (self):
        """ Tests if the module can correctly read a prior definition 
            file. """
        priors = read_priors_file ('input/simple_enzymatic.priors')
        for x in priors:
            distribution = x.get_distribution ()
            if x.name == 'k1':
                self.assertEqual (distribution.get_a (), 2.0)
                self.assertEqual (distribution.get_b (), 0.01)
            elif x.name == 'd1' or x.name == 'kcat':
                self.assertEqual (distribution.get_a (), 2.0)
                self.assertEqual (distribution.get_b (), 0.1)
            elif x.name == "Sigma":
                pass
            else:
                self.fail ()

    def test_sigma_prior (self):
        """ Tests if the module can read the experimental error 
            prior. """
        priors = read_priors_file ('input/simple_enzymatic.priors')
        sigma = priors.get_experimental_error ()
        assert (sigma > 0)


    def test_with_sbml (self):
        """ Tests if the module can define the prior of sbml parameters.
        """
        model = SBML ()
        model.load_file ("input/model1.xml")
        theta = define_sbml_params_priors (model, 'input/model1.priors')
        self.assertEqual (theta.get_size () - 1, 
                len (model.get_all_param ()))
        sigma = theta.get_experimental_error ()
        assert (sigma > 0)
