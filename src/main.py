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
        
        if which_experiment == 1:
            rand_param = RandomParameter (param, 4, .5)
        else:
            if re.search ("Km", param_original_name):
                rand_param = RandomParameter (param, 2.0, 3333.0)
            else:
                rand_param = RandomParameter (param, 1.1, 9.0)
        theta_prior.append (rand_param)
    
    if which_experiment == 1:
        sigma = RandomParameter ("experimental_sigma", 1, 1)
    else:
        sigma = RandomParameter ("experimental_sigma", 2.0, 3333.0)

    theta_prior.set_experimental_error (sigma)
    return theta_prior


which_experiment = int (argv[1])
sbml_file = argv[2]

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
    theta_priors = read_priors_file ('../input/simple_enzymatic/' + \
            'simple_enzymatic.priors')
    log_l = ml.estimate_marginal_likelihood (experiments, odes, 
            theta_priors)
    print ("log_l = " + str (log_l))
    

elif which_experiment == 1:
    sbml = SBML ()
    sbml.load_file (sbml_file)
    odes = sbml_to_odes (sbml)
    experiments = read_data_experiment_file ('../input/goodwin3.data', 
            'x1')
    theta_priors = get_theta (sbml, odes, which_experiment)
    ml = MarginalLikelihood (20000, 1000, 500, 500, 10, 10)
    log_l = ml.estimate_marginal_likelihood (experiments, odes, 
            theta_priors)
    print ("log_l = " + str (log_l))


elif which_experiment == 2:
    sbml = SBML ()
    sbml.load_file (sbml_file)
    odes = sbml_to_odes (sbml)
    t = np.linspace (0, 120, 200)
    odes.overtime_plot (['EGF', 'ERK', 'ERKPP'], t)

    experiments = []
    for i in range (1, 25):
        ex = read_data_experiment_file ('../input/Kolch/ex_' + str (i) +
            '.data', 'ERKPP')[0]
        experiments.append (ex)

    theta_priors = get_theta (sbml, odes, which_experiment)
    ml = MarginalLikelihood (50000, 1000, 10000, 2000, 10, 10)
    log_l = ml.estimate_marginal_likelihood (experiments, odes, 
            theta_priors)
    print ("log_l = " + str (log_l))


else:
    sbml = SBML ()
    sbml.load_file (sbml_file)
    odes = sbml_to_odes (sbml)
    
    # t = [0, 2, 5, 10, 20, 40, 60, 100]
    # odes.overtime_plot (['Rpp'], t)

    experiments = []
    for i in range (3):
        ex = read_data_experiment_file ('../input/bioinformatics/ex_' +
                str (i) + '.data')[0]
        experiments.append (ex)

    theta_priors = read_priors_file ('../input/bioinformatics/' + \
            'model.priors')

    ml = MarginalLikelihood (4000, 1000, 1000, 500, 20, 2)
    log_l = ml.estimate_marginal_likelihood (experiments, odes,
            theta_priors)
