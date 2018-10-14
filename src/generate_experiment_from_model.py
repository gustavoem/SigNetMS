from SBML import SBML
from ODES import ODES
from SBMLtoODES import sbml_to_odes
from Experiment import Experiment
import numpy as np


def add_noise (values):
    for i in range (len (values)):
        eps = np.random.normal (0, .1)
        if values[i] + eps > 0:
            values[i] += eps


sbml = SBML ()
sbml.load_file ('../input/simple_enzymatic.xml')
odes = sbml_to_odes (sbml)
time = np.linspace (0, 100, 6)
values = odes.evaluate_on (time)

for i in range (4):
    noised_values = {}
    for x in values:
        noised_values[x] = list (values[x])

    add_noise (noised_values["E"])
    experiment = Experiment (time, noised_values["E"], "E")
    experiment.save_to_file ('../input/simple_enzymatic_' + str (i) + 
            '.data')
