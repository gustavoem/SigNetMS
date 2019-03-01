import sys
sys.path.insert (0, '../..')

from marginal_likelihood.SBML import SBML
from marginal_likelihood.ODES import ODES
from marginal_likelihood.SBMLtoODES import sbml_to_odes
from experiment.ExperimentSet import ExperimentSet
from experiment.Experiment import Experiment
import numpy as np


def add_noise (values):
    for i in range (len (values)):
        eps = np.random.normal (0, 3)
        if values[i] + eps > 0:
            values[i] += eps

sbml = SBML ()
sbml.load_file ('final_model.sbml')
odes = sbml_to_odes (sbml)
time = [30, 60, 180, 300, 900, 1800]


# Simple experiment: run final_model simulations adding a Gaussian noise
values = odes.evaluate_exp_on ('MAPK_PP + MAPK_P', time)
experiment_set = ExperimentSet ()
for i in range (3):
    noised_values = list (values)
    add_noise (noised_values)
    experiment = Experiment (time, noised_values, 
            'MAPK_PP + MAPK_P')
    experiment_set.add (experiment)
experiment_set.save_to_file ('gauss_noise.data')


# More complex experiment: run final_model with perturbations on 
# catalytic constants
perturbation_exp = ExperimentSet ()

# First perturbation experiments
theta = odes.get_all_parameters ()
changed_param = ''
changed_param_value = 0
for p in theta:
    original_name = sbml.get_original_param_name (p)
    # kcat9 is ERK dephosphorylation catalytic constant
    if original_name == 'kcat9':
        original_value = odes.param_table[p]
        new_value = .25 * original_value
        odes.define_parameter (p, new_value)
        changed_param_value = original_value
        changed_param = p
values = odes.evaluate_exp_on ('MAPK_PP + MAPK_P', time)
for i in range (3):
    noised_values = list (values)
    add_noise (noised_values)
    experiment = Experiment (time, noised_values, 'MAPK_PP + MAPK_P')
    perturbation_exp.add (experiment)
odes.define_parameter (changed_param, changed_param_value)

# Second perturbation experiments
theta = odes.get_all_parameters ()
changed_param = ''
changed_param_value = 0
for p in theta:
    original_name = sbml.get_original_param_name (p)
    # kcat9 is MEK dephosphorylation catalytic constant
    if original_name == 'kcat5':
        original_value = odes.param_table[p]
        new_value = .25 * original_value
        odes.define_parameter (p, new_value)
        changed_param_value = original_value
        changed_param = p
        values = odes.evaluate_exp_on ('MAPK_PP + MAPK_P', time)
for i in range (3):
    noised_values = list (values)
    add_noise (noised_values)
    experiment = Experiment (time, noised_values, 'MAPK_PP + MAPK_P')
    perturbation_exp.add (experiment)

perturbation_exp.save_to_file ('perturbations.data')
