from SBML import SBML
from ODES import ODES
from SBMLtoODES import sbml_to_odes
from ExperimentSet import ExperimentSet
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
args = parser.parse_args ()
sbml_file = args.model
priors_file = args.priors
experiment_file = args.experiment

print  ("Performing marginal likelihood calculations of model: " + \
        sbml_file)


sbml = SBML ()
sbml.load_file (sbml_file)
odes = sbml_to_odes (sbml)
experiments = ExperimentSet (experiment_file)
theta_priors = define_sbml_params_priors (sbml, priors_file)

ml = MarginalLikelihood (5000, 1000, 2000, 1000, 20, 2)
log_l = ml.estimate_marginal_likelihood (experiments, odes, 
        theta_priors)
print ("log_l = " + str (log_l))
