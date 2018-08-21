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

if __name__ == '__main__':
    unittest.main ()
