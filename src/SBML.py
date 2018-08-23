import libsbml

class SBML:
    """ This class contains a representation for SBML objects. """

    def __init__ (self):
        """ Default constructor. """
        self.sbml_obj = None

    
    def load_file (self, file_name):
        """ Given an xml file construct an sbml object. """
        reader = libsbml.SBMLReader ()
        sbmldoc = reader.readSBML (file_name)
        if SBML.__check_sbml_parsing_err (sbmldoc):
            return False
        self.sbml_obj = sbmldoc
        return True
        
            
    def get_species_list (self):
        """ Returns a list with all species. """
        model = self.sbml_obj.model
        number_of_species = model.getNumSpecies ()
        
        species = []
        for i in range (number_of_species):
            species.append (model.getSpecies (i).getName ())

        return species
        return True


    def get_species_kinetic_law (self, species_name):
        """ Returns the kinetic laws of a species. """
        formula = ""
        reactions = self.__get_reactions_involving (species_name)
        for reac in reactions:
            products = [x.species for x in reac.getListOfProducts ()]
            if species_name in products:
                formula += "+ "
            else:
                formula += "- "
            formula += SBML.__reaction_rate_formula (reac)
            formula += " "
        return formula
            
    
    def get_initial_concentration (self, species_name):
        """ Returns the initial concentrations of each chemical species. 
        """
        model = self.sbml_obj.model
        species = model.getSpecies (species_name)
        initial_concentration = species.getInitialAmount ()
        return initial_concentration


    def get_name (self):
        """ Returns model name. """
        return self.sbml_obj.model.getName ()
    

    def __get_reactions_involving (self, species_name):
        """ Returns a list with all the reactions that contains a 
            species either as reactant or product. """
        model = self.sbml_obj.model
        all_reactions = model.getListOfReactions ()
        participating_reactions = []
        for reac in all_reactions:
            reacs_n_prod_ref = list (reac.getListOfReactants ()) + \
                            list (reac.getListOfProducts ())
            reacs_n_prod = [s.species for s in reacs_n_prod_ref]
            for species in reacs_n_prod:
                if species == species_name:
                    participating_reactions.append (reac)
        return participating_reactions


    @staticmethod
    def __reaction_rate_formula (reaction):
        """ Given a reaction, returns the string with the formula of its 
        reaction rate. """
        kinetic = reaction.getKineticLaw ()
        formula = kinetic.getFormula ()
        params = kinetic.getListOfParameters ()
        for param in params:
            param_name = param.getName ()
            param_value = str (param.getValue ())
            formula = formula.replace (param_name, param_value)

        return formula


    @staticmethod
    def __check_sbml_parsing_err (sbmldoc):
        """ Verifies if an SBMLReader () was succesful. """
        if sbmldoc.getNumErrors () > 0:
            sbmldoc.printErrors ()
            return True
        return False
