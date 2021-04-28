import libsbml
import re
import warnings
from model.Reaction import Reaction

class SBML:
    """ This class contains a representation for SBML objects. 
        
        This class simplifies the libsbml document to our use.

        Attributes
            sbml_obj (lisbml.SBMLDocument): a lisbml object that
                contains the information read from an sbml file.
            __parameter_values (dict): a dictionary with parameter names
                and values. The names of the parameters are simplified 
                and are different from the ones seem on the original 
                sbml document.
            __global_param (list): a list of names of global parameters 
                defined on the  sbml document.
            __local_param (dict): a dictionary in which keys are
                reaction names. The values for an entry is the
                corresponding list of reaction rate constants names.
            num_params (int): stores the number of parameters
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


    def get_copy (self):
        """ Creates a copy of this object and returns it. 
        
        Returns
            sbml_copy - an SBML object that is a copy of self.
        """
        sbml_copy = SBML  ()
        sbml_copy.sbml_obj = self.sbml_obj and self.sbml_obj.clone ()
        # pylint: disable=protected-access
        sbml_copy.__parameter_values = dict (self.__parameter_values)
        sbml_copy.__global_param = list (self.__global_param)
        sbml_copy.__local_param = dict (self.__local_param)
        sbml_copy.num_params = self.num_params
        sbml_copy.name = self.name

        return sbml_copy 

    
    def load_file (self, file_name):
        """ Loads an xml file that contains an SBML model definition.
        
            Parameters:
                file_name: a string with the file path.

            Returns
                True if parsing was successful, and False otherwise.
        """
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
        return True
        

    def write_sbmldoc_to_file (self, filename):
        """ Writes the libsbml object as an SBML file. 
        
        Parameters
            filename: a string with the name of the file to be written.
        """
        writer = libsbml.SBMLWriter ()
        writer.writeSBML (self.sbml_obj, filename)

            
    def get_species_list (self):
        """ Gets the list of species of the current model.
        
        Returns
            species: a list of strings containing species ids.
        """
        if self.sbml_obj == None:
            return []

        model = self.sbml_obj.model
        number_of_species = model.getNumSpecies ()
        
        species = []
        for i in range (number_of_species):
            species.append (model.getSpecies (i).getId ())

        return species


    def get_species_kinetic_law (self, species_name):
        """ Gets the kinetic law of a species from the current model. 
            
            Parameters
                species_name: a string with the species id.

            Returns a string wit a formula that represents the kinetic
            law that rules for the specified species.
        """
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
        """ Gets the initial concentration for a chemical species.
    
            Parameters
                species_name: the id of the species.

            Returns a number that represents the initial concentration
            of the corresponding chemical species.
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
        """ Gets the value for some model parameter.
        
            Parameters
                param: a string with the name of the parameter. Note
                that this name is most likely not the same name as used
                on the original .sbml file.

            Returns a number with the value of the corresponding
            parameters.
        """
        return self.__parameter_values[param]


    def get_all_param (self):
        """ Gets all parameters.
        
            Returns a dictionary with parameter names and values. Note
            that parameter names are likely not the same as used on the
            original .sbml file.
        """
        return self.__parameter_values


    def get_name (self):
        """ Gets model name. 
            
            Returns a string with the model name.
        """
        return self.sbml_obj.model.getName ()


    def get_original_param_name (self, param):
        """ Gets the original name of a model parameter.

            Parameters
                param: a string with a parameter name. This name is the
                same as used in self.__parameter_values, and it was
                chosen by this object using the method 
                __new_parameter ().
        
            Returns the original name of the specified parameter, that
            is, the name used in the sbml file to represent such 
            parameter. 
        """
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
        """ Gets a list with all reaction rate formulas.

            Returns
                all_formulas: a list of string containing the reaction
                    rate of all reactions.
        """
        model = self.sbml_obj.model
        reaction_list = model.getListOfReactions ()
        get_reac_f = lambda reac: reac.getKineticLaw ().getFormula ()
        return [get_reac_f (reac) for reac in reaction_list]


    def get_all_reactions (self):
        """ Gets a list with all reactions of the model. 
        
            Returns
                all_reactions: a list of Reaction objects that represent 
                    all reactions of this sbml model.
        """
        model = self.sbml_obj.model
        all_reactions = []
        for sbml_reaction in model.getListOfReactions ():
            reaction = Reaction.create_from_SBML (sbml_reaction)
            all_reactions.append (reaction)
        return all_reactions


    def set_name (self, name):
        """ Changes the name of an SBML model. 

            This method changes the attribute 'name' of the object and
            it also changes the attribute 'name' of the libsbml.Model 
            attribute.
        
            Parameters
                name: a string containing the new name.
        """
        self.name = name
        self.sbml_obj.model.setName (name)


    def add_reaction (self, reaction):
        """ Adds a new reaction to the sbml model.
        
            Parameters
                reaction: a Reaction object with the reaction to be 
                    added.

            Returns a list of string containing names of species that 
            were added to the model when adding the reaction.
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
        """ Removes a reaction from the SBML model.

            Parameters
                reaction_id: a string with the id of the reaction that
                    must be removed. If there is no reaction with id 
                    reaction_id, then the method rises a Warning and 
                    does nothing else.
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

            Returns species, a libsbml.Species object that has id 
            species_id.
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
                species: a string with the name and id of the species to 
                    be added.
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
        """ Gets all reactions involving some chemical species, as a
            product or reactant.
        
            Parameters
                species_name: a string with the id of the chemical
                species of interest.

            Returns a list with libsbml.Reaction objects, each of these
            objects representing a reaction that has the species of
            interest as a reactant or product.
        """
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
        """ Creates a new parameter

            Returns a string with the name of the new parameter.
        """
        self.num_params += 1
        return "p" + str (self.num_params)


    def __reaction_rate_formula (self, reaction):
        """ Gets the rate formula of a reaction.
        
            Parameters
                reaction: a lisbml.Reaction object representing a
                reaction.

            Returns a string with a formula of the reaction rate.
        """
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
        """ Finds the global parameters of a model.
            
            This method is used when an sbml model is loaded.
        """
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
        """ Finds the local parameters of a model.

            This method is used when an sbml model is loaded.
        """
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
