import sys
sys.path.insert (0, '../..')

from marginal_likelihood.SBML import SBML
from marginal_likelihood.ODES import ODES
from marginal_likelihood.SBMLtoODES import sbml_to_odes
from experiment.ExperimentSet import ExperimentSet
from experiment.Experiment import Experiment
import numpy as np


def add_noise (values):
    for i in range (len (values)):
        eps = np.random.normal (0, .01)
        if values[i] + eps > 0:
            values[i] += eps

sbml = SBML ()
sbml.load_file ('final_model.sbml')
odes = sbml_to_odes (sbml)

odes.print_equations ()

time = [30, 60, 180, 300, 900, 1800]
time = [5, 7, 10, 15, 20]

print (odes.evaluate_on (time))

values = odes.evaluate_exp_on ('(MAPK_PP + MAPK_P) / ' + \
        '(MAPK + MAPK_PP + MAPK_P)', time)

print ('(MAPK_PP + MAPK_P) / (MAPK + MAPK_PP + MAPK_P):')
print (values)
print ('')
print ('MAPK_PP + MAPK_P')
print (odes.evaluate_exp_on ('MAPK_PP + MAPK_P', time))


# experiment_set = ExperimentSet ()
# for i in range (3):
    # noised_values = list (values)
    # add_noise (noised_values)
    # experiment = Experiment (time, noised_values, 
            # 'MAPK_PP + MAPK_P')
    # experiment_set.add (experiment)
# experiment_set.save_to_file ('experiment.data')
