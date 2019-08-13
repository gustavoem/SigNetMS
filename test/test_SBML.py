import sys
sys.path.insert(0, '..')

import unittest
import re
from model.SBML import SBML
from model.Reaction import Reaction


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


    def test_lambdas_on_laws (self):
        """ Tests if it is possible to use functions on the kinetic 
            laws. """
        self.model = SBML ()
        self.model.load_file ("input/lambda_model.xml")
        law = self.model.get_species_kinetic_law ("R")
        assert ("EXPLE" not in law)


    def test_species_initial_concentration (self):
        """ Tests if SBML can return the initial concentration of a
            species """
        c = self.model.get_initial_concentration ("EGF")
        self.assertEqual (c, 1000)
        c = self.model.get_initial_concentration ("MEK")
        self.assertEqual (c, 3000)
    

    def test_species_initial_amount (self):
        """ Tests if SBML can return the initial concentration of a 
            species even whe it was defined as initial amount. """
        self.model = SBML ()
        self.model.load_file ("input/simple_enzymatic.xml")
        c = self.model.get_initial_concentration ("E")
        self.assertEqual (c, 10)
        c = self.model.get_initial_concentration ("S")
        self.assertEqual (c, 100)


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


    def test_get_original_parameter_name (self):
        """ Given a parameter of the our sbml model, we should be able
            to know its original name on the sbml file. """
        law = self.model.get_species_kinetic_law ("EGF")
        m = re.search (r'\- (([A-z]|_)\w*) \* EGF \* unboundEGFR', law)
        p_str = m.group (1)
        original_name = self.model.get_original_param_name (p_str)
        self.assertEqual (original_name, "k1")

        law = self.model.get_species_kinetic_law ("inactiveSos")
        m = re.search (r'\- (([A-z]|_)\w*) \* ERKPP \* inactiveSos', law)
        p_str = m.group (1)
        original_name = self.model.get_original_param_name (p_str)
        self.assertEqual (original_name, "SosKcat")


    def test_rate_when_reactant_and_product (self):
        """ If a chemical species is both reactant and product of a 
            non-reversible reaction, then the reaction shouldn't 
            contribute to the derivative of the concentration of the
            species. """
        law = self.model.get_species_kinetic_law ("ERKPP")
        self.assertFalse (re.search ("activeSos", law))
    

    def test_unnamed_parameters (self):
        """ If there's a parameter without a name, we should use its id
            as a name instead. """
        model = SBML ()
        model.load_file ("input/goodwin3.xml")
        params = model.get_all_param ()
        self.assertEqual (len (params), 4)


    def test_get_original_name_with_id (self):
        """ If there's a parameter without a name, when you call
            get_original_param_name you should get the parameter's 
            id. """
        model = SBML ()
        model.load_file ("input/goodwin3.xml")
        params = model.get_all_param ()
        names = []
        for p in params:
            original_name = model.get_original_param_name (p)
            names.append (original_name)
        
        assert ("m" in names)
        assert ("k1" in names)
        assert ("k2" in names)
        assert ("k3" in names)


    def test_add_reaction (self):
        """ Tests if it is possible to add a reaction to an SBML model. 
        """
        model = SBML ()
        model.load_file ("input/model1_bioinformatics.xml")
        
        parameters = [{"name": "kcat", "value": .5},
                      {"name": "Km", "value": 5}]
        new_reaction = Reaction ("R --dS--> Rpp", ["R"], ["Rpp"], 
                ["dS"], parameters, "kcat * dS * R / (Km + R)")
        model.add_reaction (new_reaction)
        all_formula = [f for f in model.get_all_reaction_formulas ()]
        assert (new_reaction.formula in all_formula)


if __name__ == '__main__':
    unittest.main ()
