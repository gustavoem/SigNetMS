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
    # time = np.linspace (0, 500, 100)
    # odes.overtime_plot (["E", "S", "ES", "P"], time)
    ex0 = read_data_experiment_file ('../input/' \
            'simple_enzymatic_0.data', 'E')[0]
    ex1 = read_data_experiment_file ('../input/' \
            'simple_enzymatic_1.data', 'E')[0]
    ex2 = read_data_experiment_file ('../input/' \
            'simple_enzymatic_2.data', 'E')[0]
    ex3 = read_data_experiment_file ('../input/' \
            'simple_enzymatic_3.data', 'E')[0]
    experiments = [ex0, ex1, ex2, ex3]
    ml.estimate_marginal_likelihood (experiments, sbml, odes)

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
