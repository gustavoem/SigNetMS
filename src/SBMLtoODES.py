from ODES import ODES

def sbml_to_odes (sbml):
    """ This function transforms an SBML model into an ODES model. """
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

    return odes
