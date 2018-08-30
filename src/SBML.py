import libsbml

class SBML:
    """ This class contains a representation for SBML objects. """

    def __init__ (self):
        """ Default constructor. """
        self.sbml_obj = None
        
        # Stores the values of internal parameters. 
        # This is set when calling __get_global_params and 
        # __get_local_params ()
        self.__parameter_values = {}

        # Stores internal name of global parameters.
        # The internal name of the i-th global parameter is the i-th 
        # element of this list
        self.__global_param = []

        # Stores the internal name of local parameters. 
        # __local_param is a hash of lists. The key of the hash is 
        # the name of a reaction. The internal name of the i-th local 
        # parameter of the kinetic law of this reaction is the i-th 
        # element of __local_param[reaction].
        self.__local_param = {}

        # Stores the number of parameters
        self.num_params = 0
        
                
    
    def load_file (self, file_name):
        """ Given an xml file construct an sbml object. """
        reader = libsbml.SBMLReader ()
        sbmldoc = reader.readSBML (file_name)
        if SBML.__check_sbml_parsing_err (sbmldoc):
            return False
        self.sbml_obj = sbmldoc
        self.__global_param = self.__get_global_params ()
        self.__local_param = self.__get_local_params ()
        return True
        
            
    def get_species_list (self):
        """ Returns a list with all species. """
        model = self.sbml_obj.model
        number_of_species = model.getNumSpecies ()
        
        species = []
        for i in range (number_of_species):
            species.append (model.getSpecies (i).getName ())

        return species


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
            formula += self.__reaction_rate_formula (reac)
            formula += " "
        return formula
            
    
    def get_initial_concentration (self, species_name):
        """ Returns the initial concentrations of each chemical species. 
        """
        model = self.sbml_obj.model
        species = model.getSpecies (species_name)
        initial_concentration = species.getInitialAmount ()
        return initial_concentration


    def get_param_value (self, param):
        """ Returns the value of a parameter. """
        return self.__parameter_values[param]


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

    
    def __new_parameter (self):
        self.num_params += 1
        return "p" + str (self.num_params)


    def __reaction_rate_formula (self, reaction):
        """ Given a reaction, returns the string with the formula of its 
        reaction rate. """
        kinetic = reaction.getKineticLaw ()
        formula = kinetic.getFormula ()
        params = kinetic.getListOfParameters ()
        internal_params = self.__local_param[reaction.getName ()]
        
        for i in range (len (params)):
            param_name = params[i].getName ()
            internal_param_name = internal_params[i]
            formula = formula.replace (param_name, internal_param_name)
        return formula


    def __get_global_params (self):
        model = self.sbml_obj.model
        params = model.getListOfParameters ()
        global_params = []
        for param in params:
            new_param_name = self.__new_parameter ()
            param_value = str (param.getValue ())
            global_params.append (new_param_name)
            self.__parameter_values[new_param_name] = param.getValue ()
        return global_params


    def __get_local_params (self):
        model = self.sbml_obj.model
        all_reactions = model.getListOfReactions ()
        local_params = {}
        for reac in all_reactions:
            reac_name = reac.getName ()
            local_params[reac_name] = []
            kinetic = reac.getKineticLaw ()
            params = kinetic.getListOfParameters ()
            for param in params:
                new_param_name = self.__new_parameter ()
                local_params[reac_name].append (new_param_name)
                self.__parameter_values[new_param_name] = param.getValue ()
        return local_params



    @staticmethod
    def __check_sbml_parsing_err (sbmldoc):
        """ Verifies if an SBMLReader () was succesful. """
        if sbmldoc.getNumErrors () > 0:
            sbmldoc.printErrors ()
            return True
        return False
