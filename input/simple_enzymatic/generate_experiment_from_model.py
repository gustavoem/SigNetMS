import sys
sys.path.insert (0, '../../src/')

from SBML import SBML
from ODES import ODES
from SBMLtoODES import sbml_to_odes
from ExperimentSet import ExperimentSet
from Experiment import Experiment
import numpy as np


def add_noise (values):
    for i in range (len (values)):
        eps = np.random.normal (0, .01)
        if values[i] + eps > 0:
            values[i] += eps

model_file = 'simple_enzymatic.xml'
time_array = [x for x in range (0, 101, 20)]
repetitions = 5
variable = 'E'
experiment_file = 'simple_enzymatic.data'

sbml = SBML ()
sbml.load_file (model_file)
odes = sbml_to_odes (sbml)
time = time_array
values = odes.evaluate_on (time)

experiment_set = ExperimentSet ()
for i in range (repetitions):
    noised_values = {}
    for x in values:
        noised_values[x] = list (values[x])

    add_noise (noised_values[variable])
    experiment = Experiment (time, noised_values[variable], variable)
    experiment_set.add (experiment)

experiment_set.save_to_file (experiment_file)
