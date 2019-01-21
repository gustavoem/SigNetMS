from marginal_likelihood.SBML import SBML
from marginal_likelihood.ODES import ODES
from marginal_likelihood.SBMLtoODES import sbml_to_odes
from marginal_likelihood.MarginalLikelihood import MarginalLikelihood
from marginal_likelihood.PriorsReader import read_priors_file
from marginal_likelihood.PriorsReader import define_sbml_params_priors
from experiment.ExperimentSet import ExperimentSet
import numpy as np
import re

import argparse

parser = argparse.ArgumentParser ()
parser.add_argument ("model", help="SBLM file with model definition.")
parser.add_argument ("priors", help="An XML file with the priors for" \
        + " the model parameters.")
parser.add_argument ("experiment", help="An XML file with the" \
        + " experiments observations.")
parser.add_argument ("first_sampling_iterations", help="How many" \
        + " iterations should be performed on the first step of the" \
        + " parameter sampling.")
parser.add_argument ("sigma_update_n", help="Iterations between" \
        + " each time the jumping variance is updated according to" \
        + " acceptance rate.")
parser.add_argument ("second_sampling_iterations", help="How many" \
        + " iterations should be performed on the second step of the" \
        + " parameter sampling.")
parser.add_argument ("third_sampling_iterations", help="How many" \
        + " iterations should be performed on the third step of the" \
        + " parameter sampling.")
args = parser.parse_args ()

# Problem input
sbml_file = args.model
priors_file = args.priors
experiment_file = args.experiment

# Algorithm parameters
first_step_n = int (args.first_sampling_iterations)
sigma_update_n = int (args.sigma_update_n)
second_step_n = int (args.second_sampling_iterations)
third_step_n = int (args.third_sampling_iterations)

print  ("Performing marginal likelihood calculations of model: " + \
        sbml_file)

sbml = SBML ()
sbml.load_file (sbml_file)
odes = sbml_to_odes (sbml)
experiments = ExperimentSet (experiment_file)
theta_priors = define_sbml_params_priors (sbml, priors_file)

ml = MarginalLikelihood (first_step_n, 
                         sigma_update_n, 
                         second_step_n, 
                         third_step_n, 20, 2)
log_l = ml.estimate_marginal_likelihood (experiments, odes, 
        theta_priors)
print ("log_l = " + str (log_l))
