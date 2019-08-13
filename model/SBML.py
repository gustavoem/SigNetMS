import libsbml
import re

class SBML:
    """ This class contains a representation for SBML objects. 
    """

    def __init__ (self):
        """ Default constructor. """
        self.sbml_obj = None
        
        # Stores the values of internal parameters. 
        # This is set when calling __get_global_params and 
        # __get_local_params ()
        self.__parameter_values = {}

        # Stores names of global parameters. The names of internal 
        # global parameters match the SBML name.
        self.__global_param = []

        # Stores the internal name of local parameters. 
        # __local_param is a hash of lists. The key of the hash is 
        # the name of a reaction. The internal name of the i-th local 
        # parameter of the kinetic law of this reaction is the i-th 
        # element of __local_param[reaction].
        self.__local_param = {}

        # Stores the number of parameters
        self.num_params = 0

        # Stores the model name
        self.name = ''
                
    
    def load_file (self, file_name):
        """ Given an xml file construct an sbml object. """
        reader = libsbml.SBMLReader ()
        sbmldoc = reader.readSBML (file_name)
        if SBML.__check_sbml_parsing_err (sbmldoc):
            return False

        # Converting function definitions
        converter = libsbml.SBMLFunctionDefinitionConverter ()
        converter.setDocument (sbmldoc)
        converter.convert ()

        self.sbml_obj = sbmldoc
        self.name = str (sbmldoc.model.name)
        self.__global_param = self.__get_global_params ()
        self.__local_param = self.__get_local_params ()
        return True
        
            
    def get_species_list (self):
        """ Returns a list with all species. """
        model = self.sbml_obj.model
        number_of_species = model.getNumSpecies ()
        
        species = []
        for i in range (number_of_species):
            species.append (model.getSpecies (i).getId ())

        return species


    def get_species_kinetic_law (self, species_name):
        """ Returns the kinetic laws of a species. """
        formula = ""
        reactions = self.__get_reactions_involving (species_name)
        for reac in reactions:
            products = [x.species for x in reac.getListOfProducts ()]
            reactants = [x.species for x in reac.getListOfReactants ()]
            if species_name in products and species_name in reactants:
                continue

            if species_name in products:
                formula += "+ "
            else:
                formula += "- "
            formula += self.__reaction_rate_formula (reac)
            formula += " "
        if formula == "":
            formula = "0"
        
        return SBML.__remove_compartments (formula)
            
    
    def get_initial_concentration (self, species_name):
        """ Returns the initial concentrations of each chemical species. 
        """
        model = self.sbml_obj.model
        species = model.getSpecies (species_name)
        initial_concentration = 0.0
        if species.isSetInitialConcentration ():
            initial_concentration = species.getInitialConcentration ()
        elif species.isSetInitialAmount ():
            initial_concentration = species.getInitialAmount ()
        return initial_concentration


    def get_param_value (self, param):
        """ Returns the value of a parameter. """
        return self.__parameter_values[param]


    def get_all_param (self):
        """ Return a hash with all (internal) parameters and values. """
        return self.__parameter_values


    def get_name (self):
        """ Returns model name. """
        return self.sbml_obj.model.getName ()


    def get_original_param_name (self, param):
        """ Given a parameter name (from the inside scope), returns the
            original name of this parameter inside the SBML file. """
        if param in self.__global_param:
            return param

        model_reactions = self.sbml_obj.model.getListOfReactions ()
        for reac in model_reactions:
            reac_name = reac.getId ()
            reac_internal_params = self.__local_param[reac_name]
            for i in range (len (reac_internal_params)):
                if reac_internal_params[i] == param:
                    kin_law = reac.getKineticLaw ()
                    reac_params = kin_law.getListOfParameters ()
                    original_name = reac_params[i].getId ()
                    if original_name == "":
                        original_name = reac_params[i].getId ()
                    return original_name

        return ""


    def add_reaction (self, reaction):
        """ Adds a new reaction to the sbml model.
        
        Parameters
            reaction: a Reaction object with the reaction to be added.

        Returns
            new_species: a list of string containing names of 
                species that were added to the model when adding the
                reaction.
        """
        model = self.sbml_obj.model
        created_reac = model.createReaction ()
        new_species = []
        
        for reactant in reaction.reactants:
            reactant_ref = model.getSpeciesReference (reactant)
            if reactant_ref is None:
                self.__create_new_species (reactant)
                new_species.append (reactant)
            reactant_ref = model.getSpeciesReference (reactant)
            created_reac.addReactant (reactant_ref)

        for product in reaction.products:
            product_ref = model.getSpeciesReference (product)
            if product_ref is None:
                self.__create_new_species (product)
                new_species.append (product)
            product_ref = model.getSpeciesReference (product)
            created_reac.addProduct (product_ref)

        for modifier in reaction.modifiers:
            modifier_ref = model.getSpeciesReference (modifier)
            if modifier_ref is None:
                self.__create_new_species (modifier)
                new_species.append (modifier)
            modifier_ref = model.getSpeciesReference (modifier)
            created_reac.addModifier (product_ref)
        
        created_kinetic_law = created_reac.createKineticLaw ()
        kinetic_law.setFormula = reaction.formula
        for param in reaction.params:
            created_param = kinetic_law.createParamete ()
            created_param.setId (param["name"])
            created_param.setName (param["name"])
            created_param.seValue (param["value"])
        

    def __create_new_species (self, species):
        """ Adds a new species to the model.
            
        Parameters
            species: a string with the name and id of the species to be 
                added.
        """
        model = self.sbml_obj.model
        created_reactant = model.createSpecies ()
        if created_reactant.setId (reactant):
            raise ValueError ("Could not set the id", species, 
                    "as a new species id.")
        created_reactant.setName (reactant)
        created_reactant.setInitialConcentration (0)
        created_reactant.setConstant (False)


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

    
    def __new_parameter (self):
        self.num_params += 1
        return "p" + str (self.num_params)


    def __reaction_rate_formula (self, reaction):
        """ Given a reaction, returns the string with the formula of its 
        reaction rate. """
        kinetic = reaction.getKineticLaw ()
        formula = kinetic.getFormula ()
        params = kinetic.getListOfParameters ()
        internal_params = self.__local_param[reaction.getId ()]
        
        for i in range (len (params)):
            param_name = params[i].getId ()
            if param_name is '':
                param_name = params[i].getId ()
            internal_param_name = internal_params[i]
            formula = formula.replace (param_name, internal_param_name)
        return formula


    def __get_global_params (self):
        model = self.sbml_obj.model
        params = model.getListOfParameters ()
        global_params = []
        for param in params:
            param_value = param.getValue ()
            param_name = param.getId ()
            if param_name is '':
                param_name = param.getId ()
            global_params.append (param_name)
            self.__parameter_values[param_name] = param_value
        return global_params


    def __get_local_params (self):
        model = self.sbml_obj.model
        all_reactions = model.getListOfReactions ()
        local_params = {}
        for reac in all_reactions:
            reac_name = reac.getId ()
            local_params[reac_name] = []
            kinetic = reac.getKineticLaw ()
            params = kinetic.getListOfParameters ()
            for param in params:
                new_param_name = self.__new_parameter ()
                local_params[reac_name].append (new_param_name)
                self.__parameter_values[new_param_name] = \
                        param.getValue ()
        return local_params


    @staticmethod
    def __check_sbml_parsing_err (sbmldoc):
        """ Verifies if an SBMLReader () was succesful. """
        if sbmldoc.getNumErrors () > 0:
            sbmldoc.printErrors ()
            return True
        return False


    @staticmethod
    def __remove_compartments (formula):
        """ See issue #6 on github. I'm not sure why some models use
            comparments. """
        new_formula = re.sub (r"compartment \*", "", formula)
        new_formula = re.sub (r"uVol \*", "", new_formula)
        return new_formula
