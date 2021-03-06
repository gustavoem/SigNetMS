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
        all_formula = model.get_all_reaction_formulas ()
        assert (new_reaction.formula in all_formula)
        
        # new products/reactants
        model = SBML ()
        model.load_file ("input/model1_bioinformatics.xml")
        new_reaction = Reaction ("Dummy ---> Dummy2", ["Dummy"], 
                ["Dummy2"], [], parameters, 
                "kcat * Dummy / (Km + Dummy)")
        model.add_reaction (new_reaction)
        all_formula = model.get_all_reaction_formulas ()
        assert (new_reaction.formula in all_formula)

        # new modifiers
        model = SBML ()
        model.load_file ("input/model1_bioinformatics.xml")
        new_reaction = Reaction ("R --Dummy--> Rpp", ["R"], ["Rpp"], 
                ["Dummy"], parameters, "kcat * Dummy * R / (Km + R)")
        model.add_reaction (new_reaction)
        all_formula = model.get_all_reaction_formulas ()
        assert (new_reaction.formula in all_formula)


    def test_add_reaction_with_former_unused_species (self):
        """ Tests if it we can add a new reaction using a species that
            was already defined in the model but was never present in 
            any reaction. """
        model = SBML ()
        model.load_file ("input/model_unused_species.xml")
        parameters = [{"name": "kcat", "value": .5},
                      {"name": "Km", "value": 5}]
        new_reaction = Reaction (\
                "unused --Rpp--> S", ["unused"], ["S"], ["Rpp"], \
                parameters, "kcat * Rpp * unused / (Km + unused)")
        nof_species_before = len (model.get_species_list ())
        model.add_reaction (new_reaction)
        nof_species_after = len (model.get_species_list ())
        self.assertEqual (nof_species_after, nof_species_before)


    def test_remove_reaction (self):
        """ Tests if it is possible to remove a reaction from an SBML
        model.
        """
        model = SBML ()
        model.load_file ("input/model1_bioinformatics.xml")
        model.remove_reaction ("reaction_0")
        removed_formula = "compartment * k1 * S"
        all_formula = model.get_all_reaction_formulas ()
        assert (removed_formula not in all_formula)


    def test_get_all_reaction (self):
        """ Tests if one can get all Reactions from an SBML model. """
        model = SBML ()
        model.load_file ("input/model1_bioinformatics.xml")
        gotten_reac_ids = [r.id for r in model.get_all_reactions ()]
        reac_ids = ["reaction_0", "reaction_1", "reaction_2", \
                "reaction_3"]
        for reac_id in reac_ids:
            assert reac_id in gotten_reac_ids

    def test_get_copy (self):
        """ Tests if one can copy the SBML object. """
        model = SBML ()
        copy = model.get_copy ()
        self.assertEqual (len (copy.get_species_list ()), 0)

        model.load_file ("input/model1_bioinformatics.xml")
        copy = model.get_copy ()
        self.assertEqual (len (copy.get_species_list ()), 5)
        self.assertEqual (copy.get_initial_concentration ('S'), 1)
        self.assertEqual (copy.get_initial_concentration ('dS'), 0)
        self.assertEqual (copy.get_initial_concentration ('R'), 1)
        self.assertEqual (copy.get_initial_concentration ('RS'), 0)
        self.assertEqual (copy.get_initial_concentration ('Rpp'), 0)

if __name__ == '__main__':
    unittest.main ()
