from SBML import SBML
from ODES import ODES
from SBMLtoODES import sbml_to_odes
from Experiment import Experiment
from ExperimentReader import read_experiment_file 
import MarginalLikelihood as ml

sbml = SBML ()
sbml.load_file ('../input/model1.xml')
odes = sbml_to_odes (sbml)
experiments = read_experiment_file ('../input/ERK_small.txt', 'ERK')

ml.estimate_marginal_likelihood (experiments, odes)
