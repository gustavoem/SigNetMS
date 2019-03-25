import sys
sys.path.insert (0, '../..')
import numpy as np

from marginal_likelihood.SBML import SBML
from marginal_likelihood.ODES import ODES
from marginal_likelihood.SBMLtoODES import sbml_to_odes

for i in range (1, 5):
    model = 'model' + str (i)
    sbml = SBML ()
    sbml.load_file (model + '.xml')
    odes = sbml_to_odes (sbml)
    time = np.linspace (0, 1000, 7)
    odes.overtime_plot (['B', 'A', 'AB', 'C'], time, filename=model)

