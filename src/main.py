from SBML import SBML
from ODES import ODES
from SBMLtoODES import sbml_to_odes
from Experiment import Experiment
from ExperimentReader import read_txt_experiment_file 
from ExperimentReader import read_data_experiment_file 
from RandomParameterList import RandomParameterList
from RandomParameter import RandomParameter
from MarginalLikelihood import MarginalLikelihood
import numpy as np
import re
from sys import argv

def get_theta (sbml, model, which_experiment):
    """ Given a model, construct a list containing all parameters of 
        the model as random variables. """
    theta_prior = RandomParameterList ()
    params = model.get_all_parameters ()
    
    for param in params:
        param_original_name = sbml.get_original_param_name (param)
        
        if which_experiment == 0:
            if param_original_name == "k1":
                rand_param = RandomParameter (param, 2.0, 0.01)
            if param_original_name == "d1" or param_original_name == "kcat":
                rand_param = RandomParameter (param, 2.0, 0.1)
        elif which_experiment == 1:
            rand_param = RandomParameter (param, 4, .5)
        else:
            if re.search ("Km", param_original_name):
                rand_param = RandomParameter (param, 2.0, 3333.0)
            else:
                rand_param = RandomParameter (param, 1.1, 9.0)
        theta_prior.append (rand_param)
    
    if which_experiment == 0:
        sigma = RandomParameter ("experimental_sigma", 2.0, 2.6)
    elif which_experiment == 1:
        sigma = RandomParameter ("experimental_sigma", 1, 1)
    else:
        sigma = RandomParameter ("experimental_sigma", 2.0, 3333.0)

    theta_prior.set_experimental_error (sigma)
    return theta_prior


which_experiment = int (argv[1])

if which_experiment == 0:
    sbml = SBML ()
    sbml.load_file ('../input/simple_enzymatic/simple_enzymatic.xml')
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
    ml = MarginalLikelihood (1000, 500, 500, 500, 10, 10)
    theta_priors = get_theta (sbml, odes, which_experiment)
    log_l = ml.estimate_marginal_likelihood (experiments, odes, 
            theta_priors)
    print ("log_l = " + str (log_l))
    
elif which_experiment == 1:
    sbml = SBML ()
    sbml.load_file ('../input/goodwin3.xml')
    odes = sbml_to_odes (sbml)
    experiments = read_data_experiment_file ('../input/goodwin3.data', 
            'x1')
    theta_priors = get_theta (sbml, odes, which_experiment)
    ml = MarginalLikelihood (20000, 1000, 500, 500, 10, 10)
    log_l = ml.estimate_marginal_likelihood (experiments, odes, 
            theta_priors)
    print ("log_l = " + str (log_l))

else:
    sbml = SBML ()
    sbml.load_file ('../input/Kolch/model2.xml')
    odes = sbml_to_odes (sbml)

    experiments = []
    for i in range (1, 25):
        ex = read_data_experiment_file ('../input/Kolch/ex_' + str (i) +
            '.data', 'ERKPP',)[0]
        experiments.append (ex)
    theta_priors = get_theta (sbml, odes, which_experiment)
    ml = MarginalLikelihood (50000, 1000, 10000, 2000, 10, 10)
    log_l = ml.estimate_marginal_likelihood (experiments, odes, 
            theta_priors)
    print ("log_l = " + str (log_l))
