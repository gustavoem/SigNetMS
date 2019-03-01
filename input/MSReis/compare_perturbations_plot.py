import sys
sys.path.insert (0, '../..')

from marginal_likelihood.SBML import SBML
from marginal_likelihood.ODES import ODES
from marginal_likelihood.SBMLtoODES import sbml_to_odes


sbml = SBML ()
sbml.load_file ('final_model.sbml')
odes = sbml_to_odes (sbml)
time = [0, 30, 60, 180, 300, 900, 1800]
odes.overtime_plot (['MAPK_PP + MAPK_P'], time, filename='unperturbed')

theta = odes.get_all_parameters ()
changed_param = ''
changed_param_value = 0
for p in theta:
    original_name = sbml.get_original_param_name (p)
    # kcat9 is ERK dephosphorylation catalytic constant
    if original_name == 'kcat9':
        original_value = odes.param_table[p]
        new_value = .1 * original_value
        odes.define_parameter (p, new_value)
        changed_param_value = original_value
        changed_param = p
odes.overtime_plot (['MAPK_PP + MAPK_P'], time, 
        filename='perturbed_erk_dephosphorylation')
odes.define_parameter (changed_param, changed_param_value)


theta = odes.get_all_parameters ()
changed_param = ''
changed_param_value = 0
for p in theta:
    original_name = sbml.get_original_param_name (p)
    # kcat5 is MEK dephosphorylation catalytic constant
    if original_name == 'kcat5':
        original_value = odes.param_table[p]
        new_value = .1 * original_value
        odes.define_parameter (p, new_value)
        changed_param_value = original_value
        changed_param = p
odes.overtime_plot (['MAPK_PP + MAPK_P'], time, 
        filename='perturbed_mek_dephosphorylation')
odes.define_parameter (changed_param, changed_param_value)
