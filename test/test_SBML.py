import sys
sys.path.insert(0, '../src/')

import unittest
import re
from SBML import SBML


class TestSBMLMethods (unittest.TestCase):

    def setUp (self):
        self.model = SBML ()
        self.model.load_file ("input/model1.xml")


    def test_SBML_parsing (self):
        """ Tests if SBML model can store the model name. """
        self.assertEqual ("Model 1: One branch", self.model.get_name ())


    def test_species_list (self):
        """ Tests if SBML stores the list of species in the model. """
        species_list = self.model.get_species_list ()
        self.assertEqual (27, len (species_list))


    def test_species_kinetic_law (self):
        """ Tests if SBML can return the equation of rate change of a 
            species. """
        law = self.model.get_species_kinetic_law ("unboundEGFR")
        m1 = re.search (r'\- (([A-z]|_)\w*) \* EGF', law)
        m2 = re.search (r'\+ (([A-z]|_)\w*) \* boundEGF', law)
        self.assertGreater (len (m1.groups ()), 0)
        self.assertGreater (len (m2.groups ()), 0)
        p1_str = m1.group (1)
        p2_str = m2.group (1)
        p1 = self.model.get_param_value (p1_str)
        p2 = self.model.get_param_value (p2_str)
        #eq = '- 1.6665105 * EGF * unboundEGFR + 687.15641 * boundEGFR'
        self.assertEqual (p1, 1.6665105)
        self.assertEqual (p2, 687.15641)


    def test_species_initial_concentration (self):
        """ Tests if SBML can return the initial concentration of a
            species """
        c = self.model.get_initial_concentration ("EGF")
        self.assertEqual (c, 1000)
        c = self.model.get_initial_concentration ("MEK")
        self.assertEqual (c, 3000)


    def test_parameter_consistency (self):
        """ Tests if parameters from reactions are consistently used on
            the equations of rate change """
        law = self.model.get_species_kinetic_law ("unboundEGFR")
        m1 = re.search (r'\- (([A-z]|_)\w*) \* EGF', law)
        p1_str = m1.group (1)
        
        law2 = self.model.get_species_kinetic_law ("boundEGFR")
        m2 = re.search (r'\+ (([A-z]|_)\w*) \* EGF', law)
        p2_str = m1.group (1)

        self.assertEqual (p1_str, p2_str)
        self.assertEqual (self.model.get_param_value (p1_str), 
                self.model.get_param_value (p2_str))


if __name__ == '__main__':
    unittest.main ()
