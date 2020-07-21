import sys
sys.path.insert (0, '../..')

from model.SBML import SBML
from model.ODES import ODES
from model.SBMLtoODES import sbml_to_odes
from experiment.ExperimentSet import ExperimentSet
from experiment.Experiment import Experiment
import numpy as np


def add_noise (values):
    for i in range (len (values)):
        eps = np.random.normal (0, .01)
        if values[i] + eps > 0:
            values[i] += eps

sbml = SBML ()
sbml.load_file ('Correct-model.sbml')
odes = sbml_to_odes (sbml)
time = [0, 30, 60, 90, 120, 150, 180, 210, 240]
values = odes.evaluate_on (time)

experiment_set = ExperimentSet ()
for i in range (3):
    noised_values = {}
    for x in values:
        noised_values[x] = list (values[x])

    add_noise (noised_values["RasGTP"])
    experiment = Experiment (time[1:], noised_values["RasGTP"][1:], 
            "RasGTP")
    experiment_set.add (experiment)

experiment_set.save_to_file ('experiment.data')
