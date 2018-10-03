from SBML import SBML
from ODES import ODES
from SBMLtoODES import sbml_to_odes
from Experiment import Experiment
from ExperimentReader import read_txt_experiment_file 
from ExperimentReader import read_data_experiment_file 
import MarginalLikelihood as ml
import numpy as np

which_experiment = 0

if which_experiment == 0:
    sbml = SBML ()
    sbml.load_file ('../input/simple_enzymatic.xml')
    odes = sbml_to_odes (sbml)
    time = np.linspace (0, 20, 11)
    odes.overtime_plot (["E", "S", "ES", "P"], time)

elif which_experiment == 1:
    sbml = SBML ()
    sbml.load_file ('../input/goodwin3.xml')
    odes = sbml_to_odes (sbml)
    experiments = read_data_experiment_file ('../input/goodwin3.data', 
            'x1')
    ml.estimate_marginal_likelihood (experiments, sbml, odes)

else:
    sbml = SBML ()
    sbml.load_file ('../input/model1.xml')
    odes = sbml_to_odes (sbml)
    experiments = read_txt_experiment_file ('../input/ERK_data.txt', 
            'ERK')
    ml.estimate_marginal_likelihood (experiments, sbml, odes)
