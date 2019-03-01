import sys
sys.path.insert (0, '../..')

from marginal_likelihood.SBML import SBML
from marginal_likelihood.ODES import ODES
from marginal_likelihood.SBMLtoODES import sbml_to_odes


sbml = SBML ()
sbml.load_file ('model2.xml')
odes = sbml_to_odes (sbml)
time = [0, 2, 5, 10, 20, 40]
time = [2, 5, 10]
# time = [0, 120, 300, 600] # in seconds

initial_state = {'ERKPP': 903, 'ERK': 9097}
odes.overtime_plot (['ERKPP/100'], time, filename='sum', 
        initial_state_map=initial_state)
