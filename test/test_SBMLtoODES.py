import sys
sys.path.insert (0, '../src/')

import unittest
from SBML import SBML
from ODES import ODES
from SBMLtoODES import sbml_to_odes 

class TestSBMLtoODES (unittest.TestCase):
    
    def test_sbml_to_odes (self):
        """ Tests if you can create an ODES object from an SBML model. 
        """
        model = SBML ()
        model.load_file ("input/model1.xml")
        odes = sbml_to_odes (model)
        y = odes.evaluate_on ([0])
        self.assertListEqual (y["EGF"], [1000])
        self.assertListEqual (y["ERK"], [10000])

