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
sbml.load_file ('initial_model.sbml')
odes = sbml_to_odes (sbml)
time = [0.5, 1, 3, 5, 15, 30]
values = odes.evaluate_on (time)

experiment_set = ExperimentSet ()
for i in range (3):
    noised_values = {}
    for x in values:
        noised_values[x] = list (values[x])

    add_noise (noised_values["MAPK_PP"])
    experiment = Experiment (time[1:], noised_values["MAPK_PP"][1:], 
            "MAPK_PP")
    experiment_set.add (experiment)

experiment_set.save_to_file ('experiment.data')
