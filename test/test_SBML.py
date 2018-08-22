import sys
sys.path.insert(0, '../src/')

import unittest
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
        """ Tests if SBML can return the kinetic law of a species. """
        law = self.model.get_species_kinetic_law ("unboundEGFR")
        eq = '- 1.6665105 * EGF * unboundEGFR + 687.15641 * boundEGFR'
        self.assertEqual (law, eq)


    def test_species_initial_concentration (self):
        """ Tests if SBML can return the initial concentration of a
            species """
        c = self.model.get_initial_concentration ("EGF")
        self.assertEqual (c, 1000)
        c = self.model.get_initial_concentration ("MEK")
        self.assertEqual (c, 3000)



if __name__ == '__main__':
    unittest.main ()
