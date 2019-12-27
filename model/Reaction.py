class Reaction:
    """ This class represents reactions that can be addded to SBML
        models.

        Attributes:
            id (string): The identification of the reaction.
            reactants (list): A list of strings containing all the
                reactants.
            products (list): A list of strings containing all the 
                products.
            modifiers (list): A list of strings containing all the 
                modifiers.
            parameters (list): A list of dictionaries containing the 
                keys: "name" and "value".
            formula (string): The formula of the rate of the reaction.
    """

    def __init__ (self, reac_id, reactants, products, modifiers, \
            parameters, formula):
        """ The constructor of Reaction.

            Parameters
                id: a string with the identification of the reaction.
                reactants: a list of strings containing all the
                    reactants.
                products: a list of strings containing all the products.
                modifiers: a list of strings containing all the
                    modifiers.
                parameters: a list of dictionaries containing the keys: 
                    "name" and "value".
                formula: The formula of the rate of the reaction.
        """
        self.id = reac_id
        self.reactants = reactants
        self.products = products
        self.modifiers = modifiers
        self.parameters = parameters
        self.formula = formula

    
    @staticmethod
    def create_from_SBML (sbml_reaction):
        """Creates a Reaction object given an libsbml.Reaction object.
        
        Parameters
            sbml_reaction: an SBML Reaction object.

        Returns
            reaction: a Reaction object that is equivalent to 
                sbml_reaction.
        """
        reac_id = sbml_reaction.getId ()
        reactants = [s.species \
                for s in sbml_reaction.getListOfReactants()]
        products = [s.species \
                for s in sbml_reaction.getListOfProducts()]
        modifiers = [s.species \
                for s in sbml_reaction.getListOfModifiers()]
        kinetic_law = sbml_reaction.getKineticLaw ()
        formula = kinetic_law.getFormula ()
        parameters = []
        for p in kinetic_law.getListOfParameters ():
            id_name = p.getName ()
            value = p.getValue ()
            param = {"name": id_name, "value": value}
            parameters.append (param)
        reaction = Reaction (reac_id, reactants, products, modifiers,\
                parameters, formula)
        return reaction
