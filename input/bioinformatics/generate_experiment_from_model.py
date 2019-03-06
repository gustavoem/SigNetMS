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
sbml.load_file ('model1.xml')
odes = sbml_to_odes (sbml)
time = [0, 2, 5, 10, 20, 40, 60, 100]
values = odes.evaluate_on (time)

experiment_set = ExperimentSet ()
for i in range (3):
    noised_values = {}
    for x in values:
        noised_values[x] = list (values[x])

    add_noise (noised_values["Rpp"])
    experiment = Experiment (time[1:], noised_values["Rpp"][1:], 
            "Rpp")
    experiment_set.add (experiment)

experiment_set.save_to_file ('experiment.data')
