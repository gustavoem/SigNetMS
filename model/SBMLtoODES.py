from model.ODES import ODES

def sbml_to_odes (sbml):
    """ Creates an ODE model given an SBML object.
    
        Parameters:
            sbml: a model.SBML object.

        Returns an ODES object that is a system of differential
        equations that rules the concentration changes of the chemical
        species specified on the sbml model.
    """
    odes = ODES ()
    variables = sbml.get_species_list ()
    for var in variables:
        formula = sbml.get_species_kinetic_law (var)
        initial_val = sbml.get_initial_concentration (var)
        odes.add_equation (var, formula)
        odes.define_initial_value (var, initial_val)

    params = sbml.get_all_param ()
    for param in params:
        odes.define_parameter (param, params[param])
    odes.name = sbml.name
    return odes
