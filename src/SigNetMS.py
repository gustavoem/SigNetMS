from SBML import SBML
from ODES import ODES
from SBMLtoODES import sbml_to_odes
from Experiment import Experiment
from ExperimentReader import read_txt_experiment_file 
from ExperimentReader import read_data_experiment_file 
from RandomParameterList import RandomParameterList
from RandomParameter import RandomParameter
from MarginalLikelihood import MarginalLikelihood
from PriorsReader import read_priors_file
from PriorsReader import define_sbml_params_priors
import numpy as np
import re

import argparse

parser = argparse.ArgumentParser ()
parser.add_argument ("model", help="SBLM file with model definition.")
parser.add_argument ("priors", help="An XML file with the priors for " \
        + "the model parameters.")
parser.add_argument ("experiment", help="An XML file with the" \
        + "experiments observations.")
parser.add_argument ("which_experiment", help="we're going to remove this.")
args = parser.parse_args ()
sbml_file = args.model
priors_file = args.priors
experiment = args.experiment
which_experiment = int (args.which_experiment)

print  ("Performing marginal likelihood calculations of model: " + \
        sbml_file)

if which_experiment == 0:
    sbml = SBML ()
    sbml.load_file (sbml_file)
    odes = sbml_to_odes (sbml)
    
    ex0 = read_data_experiment_file ('../input/simple_enzymatic/' + \
            'simple_enzymatic_0.data')[0]
    ex1 = read_data_experiment_file ('../input/simple_enzymatic/' + \
            'simple_enzymatic_1.data')[0]
    ex2 = read_data_experiment_file ('../input/simple_enzymatic/' + \
            'simple_enzymatic_2.data')[0]
    ex3 = read_data_experiment_file ('../input/simple_enzymatic/' + \
            'simple_enzymatic_3.data')[0]
    experiments = [ex0, ex1, ex2, ex3]
    ml = MarginalLikelihood (20000, 1000, 5000, 1000, 20, 2)
    theta_priors = define_sbml_params_priors (sbml, priors_file)
    log_l = ml.estimate_marginal_likelihood (experiments, odes, 
            theta_priors)
    print ("log_l = " + str (log_l))
    

elif which_experiment == 1:
    sbml = SBML ()
    sbml.load_file (sbml_file)
    odes = sbml_to_odes (sbml)
    t = np.linspace (0, 120, 200)
    odes.overtime_plot (['EGF', 'ERK', 'ERKPP'], t)

    experiments = []
    for i in range (1, 25):
        ex = read_data_experiment_file ('../input/Kolch/ex_' + str (i) +
            '.data')[0]
        experiments.append (ex)

    theta_priors = define_sbml_params_priors (sbml, priors_file)
    ml = MarginalLikelihood (50000, 1000, 10000, 2000, 10, 10)
    log_l = ml.estimate_marginal_likelihood (experiments, odes, 
            theta_priors)
    print ("log_l = " + str (log_l))


else:
    sbml = SBML ()
    sbml.load_file (sbml_file)
    odes = sbml_to_odes (sbml)
    
    experiments = []
    for i in range (3):
        ex = read_data_experiment_file ('../input/bioinformatics/ex_' +
                str (i) + '.data')[0]
        experiments.append (ex)

    theta_priors = define_sbml_params_priors (sbml, priors_file)
    ml = MarginalLikelihood (4000, 1000, 1000, 500, 20, 2)
    log_l = ml.estimate_marginal_likelihood (experiments, odes,
            theta_priors)
