import libsbml
import re
import warnings
from model.Reaction import Reaction

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

        # Stores the last loaded file
        self.__loaded_file = None


    def get_copy (self):
        """ Creates a copy of this object and returns it. 
        
        Returns
            sbml_copy - an SBML object that is a copy of self.
        """
        sbml_copy = SBML  ()
        if self.__loaded_file:
            sbml_copy.load_file (self.__loaded_file) 
        return sbml_copy 

    
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

        self.num_params = 0
        self.__parameter_values = {}
        self.sbml_obj = sbmldoc
        self.name = str (sbmldoc.model.name)
        self.__global_param = self.__get_global_params ()
        self.__local_param = self.__get_local_params ()
        self.__loaded_file = file_name
        return True
        

    def write_sbmldoc_to_file (self, filename):
        """ Writes the libsbml object as an SBML file. 
        
        Parameters
            filename: a string with the name of the file to be written.
        """
        writer = libsbml.SBMLWriter ()
        writer.writeSBML (self.sbml_obj, filename)

            
    def get_species_list (self):
        """ Returns a list with all species. """
        if self.sbml_obj == None:
            return []

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


    def get_all_reaction_formulas (self):
        """ Returns a list with all reaction rate formulas.

        Returns
            all_formulas: a list of string containing the reaction rate
                of all reactions.
        """
        model = self.sbml_obj.model
        reaction_list = model.getListOfReactions ()
        get_reac_f = lambda reac: reac.getKineticLaw ().getFormula ()
        return [get_reac_f (reac) for reac in reaction_list]


    def get_all_reactions (self):
        """ Returns a list with all reactions of the model. 
        
        Returns
            all_reactions: a list of Reaction objects that represent all
                reactions of this sbml model.
        """
        model = self.sbml_obj.model
        all_reactions = []
        for sbml_reaction in model.getListOfReactions ():
            reaction = Reaction.create_from_SBML (sbml_reaction)
            all_reactions.append (reaction)
        return all_reactions


    def set_name (self, name):
        """ Changes the name of an SBML model. 

        This method changes the attribute 'name' of the object and it
        also changes the attribute 'name' of the libsbml.Model 
        attribute.
        
        Parameters
            name: a string containing the new name.
        """
        self.name = name
        self.sbml_obj.model.setName (name)


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
        created_reac.setIdAttribute (reaction.id)
        new_species = []

        for reactant in reaction.reactants:
            model_reactant = self.__get_species (reactant)
            created_reac.addReactant (model_reactant)

        for product in reaction.products:
            model_product = self.__get_species (product)
            created_reac.addProduct (model_product)

        for modifier in reaction.modifiers:
            created_modifier = created_reac.createModifier ()
            if model.getSpeciesReference (modifier) is None:
                self.__create_new_species (modifier)
                new_species.append (modifier)
            havnt_added_mod = created_modifier.setSpecies (modifier)
            if havnt_added_mod:
                raise ValueError ("Could not set", modifier, \
                        "as a modifier")
        
        created_kinetic_law = created_reac.createKineticLaw ()
        created_kinetic_law.setFormula (reaction.formula)
        for param in reaction.parameters:
            created_param = created_kinetic_law.createParameter ()
            created_param.setIdAttribute (param["name"])
            created_param.setName (param["name"])
            created_param.setValue (param["value"])

        
    def remove_reaction (self, reaction_id):
        """
        Removes a reaction from the SBML model.

        Parameters
            reaction_id: a string with the id of the reaction that must
                be removed. If there is no reaction with id reaction_id,
                then the method rises a Warning and does nothing else.
        """
        model = self.sbml_obj.model
        reaction = model.removeReaction (reaction_id)
        if not reaction:
            warnings.warn ("Could not find reaction with id" \
                    + str (reaction_id))


    def __get_species (self, species_id):
        """ Returns species (and creates if needed) with id species_id.

        Parameters
            species_id: a string with the id of the species.

        Returns
            species: a libsbml.Species object that has id species_id.
        """
        model = self.sbml_obj.model
        species = model.getSpecies (species_id)
        if species is None:
            self.__create_new_species (species_id)
        species = model.getSpecies (species_id)
        return species


    def __create_new_species (self, species):
        """ Adds a new species to the model.
            
        Parameters
            species: a string with the name and id of the species to be 
                added.
        """
        model = self.sbml_obj.model
        created_species = model.createSpecies ()
        if created_species.setIdAttribute (species):
            raise ValueError ("Could not set the id", species, 
                    "as a new species id.")
        created_species.setName (species)
        created_species.setInitialConcentration (0)
        created_species.setConstant (False)
        # We assume there is only one compartment
        compartments = model.getListOfCompartments ()
        compartment_id = compartments[-1].getIdAttribute ()
        created_species.setCompartment (compartment_id)


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
