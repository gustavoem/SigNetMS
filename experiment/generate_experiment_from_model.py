import sys
sys.path.insert (0, '..')

from model.SBML import SBML
from model.SBMLtoODES import sbml_to_odes
from .Experiment import Experiment
from .ExperimentSet import ExperimentSet
import numpy as np


def add_noise (vals):
    for j in range (len (vals)):
        eps = np.random.normal (0, .01)
        if vals[j] + eps > 0:
            vals[j] += eps


sbml = SBML ()
sbml.load_file ('../input/bioinformatics/model1.xml')
odes = sbml_to_odes (sbml)
time = [0, 2, 5, 10, 20, 40, 60, 100]
values = odes.evaluate_on (time)

exp_set = ExperimentSet ()
for i in range (3):
    noised_values = {}
    for x in values:
        noised_values[x] = list (values[x])

    add_noise (noised_values["Rpp"])
    experiment = Experiment (time, noised_values["Rpp"], "Rpp")
    exp_set.add (experiment)

exp_set.save_to_file ('../input/bioinformatics/experiment.data')
