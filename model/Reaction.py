class Reaction:
    """ This class represents reactions that can be addded to SBML models.

    Attributes:
        name (string): The identification of the reaction.
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

    def __init__ (self, name, reactants, products, modifiers, \
            parameters, formula):
        """The constructor of Reaction.

        Parameters
            name (string): The identification of the reaction.
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
        self.name = name
        self.reactants = reactants
        self.products = products
        self.modifiers = modifiers
        self.parameters = parameters
        self.formula = formula
