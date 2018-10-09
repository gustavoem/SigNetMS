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
time = np.linspace (0, 200, 201)
values = odes.evaluate_on (time)
add_noise (values["E"])
experiment = Experiment (time, values["E"], "E")
experiment.save_to_file ('../input/simple_enzymatic.data')
