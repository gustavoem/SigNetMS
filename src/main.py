from SBML import SBML
from ODES import ODES
from SBMLtoODES import sbml_to_odes
from Experiment import Experiment
from ExperimentReader import read_txt_experiment_file 
from ExperimentReader import read_data_experiment_file 
from MarginalLikelihood import MarginalLikelihood
import numpy as np

which_experiment = 1 

if which_experiment == 0:
    sbml = SBML ()
    sbml.load_file ('../input/simple_enzymatic/simple_enzymatic.xml')
    odes = sbml_to_odes (sbml)
    # time = np.linspace (0, 500, 100)
    # odes.overtime_plot (["E", "S", "ES", "P"], time)
    ex0 = read_data_experiment_file ('../input/simple_enzymatic/' + \
            'simple_enzymatic_0.data', 'E')[0]
    ex1 = read_data_experiment_file ('../input/simple_enzymatic/' + \
            'simple_enzymatic_1.data', 'E')[0]
    ex2 = read_data_experiment_file ('../input/simple_enzymatic/' + \
            'simple_enzymatic_2.data', 'E')[0]
    ex3 = read_data_experiment_file ('../input/simple_enzymatic/' + \
            'simple_enzymatic_3.data', 'E')[0]
    experiments = [ex0, ex1, ex2, ex3]
    ml = MarginalLikelihood (5000, 2000, 1000, 10, 10)
    log_l = ml.estimate_marginal_likelihood (experiments, sbml, odes)
    print ("log_l = " + str (log_l))
    

elif which_experiment == 1:
    sbml = SBML ()
    sbml.load_file ('../input/goodwin3.xml')
    odes = sbml_to_odes (sbml)
    experiments = read_data_experiment_file ('../input/goodwin3.data', 
            'x1')
    ml = MarginalLikelihood (20000, 500, 500, 10, 20)
    ml.estimate_marginal_likelihood (experiments, sbml, odes)

else:
    sbml = SBML ()
    sbml.load_file ('../input/Kolch/model2.xml')
    odes = sbml_to_odes (sbml)
    # time = np.linspace (0, 5000, 100)
    # odes.overtime_plot (["ERK", "ERKPP"], time)

    experiments = []
    for i in range (1, 25):
        ex = read_data_experiment_file ('../input/Kolch/ex_' + str (i) +
            '.data', 'ERKPP',)[0]
        experiments.append (ex)
    ml = MarginalLikelihood (20000, 200, 100, 10, 10)
    ml.estimate_marginal_likelihood (experiments, sbml, odes)
