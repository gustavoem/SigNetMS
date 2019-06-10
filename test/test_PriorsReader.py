import sys
sys.path.insert (0, '..')

import unittest
from model.PriorsReader import read_priors_file
from model.PriorsReader import define_sbml_params_priors 
from model.SBML import SBML
import numpy as np

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
    
    
    def test_read_lognormal (self):
        """ Tests if the module can correctly read a prior definition 
            file. """
        priors = read_priors_file ('input/lognormal.priors')
        for x in priors:
            distribution = x.get_distribution ()
            if x.name == 'k1':
                expected_mean = np.exp (3)
                self.assertEqual (distribution.mean (), expected_mean)
            elif x.name == 'd1':
                expected_mean = np.exp (10)
                self.assertEqual (distribution.mean (), expected_mean)
            elif x.name == "Sigma" or x.name == "Noise":
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


    def test_definition_of_sbml_parameters (self):
        """ Tests if the priors have the correct parameter names of the
        sbml model. """
        model = SBML ()
        model.load_file ("input/model1.xml")
        theta = define_sbml_params_priors (model, 'input/model1.priors')
        sbml_parameters = model.get_all_param ()
        for t in theta.get_model_parameters ():
            assert (t.name in sbml_parameters)
        

    def test_priors_fully_defined (self):
        """ Tests if priors file define priors for all parameters of 
            a model. """
        model = SBML ()
        model.load_file ('input/model1.xml')
        self.assertRaises (ValueError, define_sbml_params_priors, model, 
                'input/incomplete.priors')

    def test_no_spurious_parameter (self):
        """ Tests if warning is showed whenever the priors file define
            some parameter that is not present on the sbml parameters.
        """
        model = SBML ()
        model.load_file ('input/model1.xml')
        define_sbml_params_priors (model, 'input/overcomplete.priors')
        self.assertWarns (Warning, define_sbml_params_priors, model, \
                'input/overcomplete.priors')
