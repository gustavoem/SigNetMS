from SBML import SBML
from ODES import ODES
from SBMLtoODES import sbml_to_odes
from Experiment import Experiment
from ExperimentReader import read_experiment_file 
import MarginalLikelihood as ml
import numpy as np

sbml = SBML ()
sbml.load_file ('../input/model1.xml')
# sbml.load_file ('../input/goodwin3.xml')
odes = sbml_to_odes (sbml)

params = odes.get_all_parameters ()
for p in params:
    print (p)

experiments = read_experiment_file ('../input/ERK_small.txt', 'ERK')
# odes.overtime_plot ("ERK", experiments[0].times)
# odes.overtime_plot ("x1", np.linspace (0, 1200, num=100))
# ml.estimate_marginal_likelihood (experiments, sbml, odes)
