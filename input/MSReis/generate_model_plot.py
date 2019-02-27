import sys
sys.path.insert (0, '../..')

from marginal_likelihood.SBML import SBML
from marginal_likelihood.ODES import ODES
from marginal_likelihood.SBMLtoODES import sbml_to_odes


sbml = SBML ()
sbml.load_file ('final_model.sbml')
odes = sbml_to_odes (sbml)
# time = [30, 60, 180, 300, 900, 1800]
time = [0, 30, 60, 180, 300, 900, 1800]
odes.overtime_plot (['MAPK_PP + MAPK_P'], time, filename='sum')
# odes.overtime_plot (['(MAPK_PP + MAPK_P) / (MAPK_PP + MAPK_P + MAPK)'], time, filename='ratio')

