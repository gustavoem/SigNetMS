from SBML import SBML
from ODES import ODES
from SBMLtoODES import sbml_to_odes
from Experiment import Experiment
import numpy as np


def add_noise (values):
    for i in range (len (values)):
        eps = np.random.normal (0, .01)
        if values[i] + eps > 0:
            values[i] += eps


sbml = SBML ()
sbml.load_file ('../input/bioinformatics/model1.xml')
odes = sbml_to_odes (sbml)
time = [0, 2, 5, 10, 20, 40, 60, 100]
values = odes.evaluate_on (time)

for i in range (3):
    noised_values = {}
    for x in values:
        noised_values[x] = list (values[x])

    add_noise (noised_values["Rpp"])
    experiment = Experiment (time, noised_values["Rpp"], "Rpp")
    experiment.save_to_file ('../input/bioinformatics/ex_' + str (i) + \
            '.data')
