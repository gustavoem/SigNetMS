import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import numpy as np
import re

sys.path.insert (0, '..')
from model.SBML import SBML
from model.SBMLtoODES import sbml_to_odes
from experiment.ExperimentSet import ExperimentSet


if len (sys.argv) != 4:
    print ("Usage: " + sys.argv[0] + " trace_file model_file " + \
            "experiment_file")
    exit()

trace_file = sys.argv[1] 
model_file = sys.argv[2]
experiment_file = sys.argv[3]

# define model
sbml = SBML ()
sbml.load_file (model_file)
ode = sbml_to_odes (sbml)
ode.print_equations ()

# get param names
param_names = []
for param in sbml.get_all_param ():
    param_names.append (param)

# experiment info
exp_set = ExperimentSet (experiment_file)
experiment_times = exp_set[-1].times
experiment_measure = exp_set[-1].measure_expression
experiment_observations = []
for exp in exp_set:
    obs = exp.values
    experiment_observations.append (obs)

parameter_regex = re.compile (r"Current theta: \[(.*)\]")

f = open (trace_file)
theta = None

step = 0
for line in f:
    m = parameter_regex.match(line)
    if m:
        theta_str = m.group (1).split (',')
        theta = np.array ([float (v) for v in theta_str])
        
    if 'Accepted' in line:
        for idx, param in enumerate (param_names):
            ode.define_parameter (param, theta[idx])
        simulation = ode.evaluate_exp_on (experiment_measure,
                experiment_times)
        
        figname = 'simulation_step_' + str (step) + '.png'
        plot_title = 'Simulation for model ' + sbml.name \
                + ', iteration ' + str (step)
        step += 1

        fig, ax =  plt.subplots ()
        # plot experiments
        i = 1
        for obs in experiment_observations:
            label = 'Experimental observation #' + str (i)
            ax.plot(experiment_times, obs, label=label)
            i += 1

        print(theta, simulation)
        ax.plot (experiment_times, simulation, label='Simulated observation')

        plt.ylabel ('$[' + experiment_measure + ']$')
        plt.xlabel ('Time (s)')
        ax.legend ()
        fig.savefig (figname, transparent=True)
        plt.clf ()
        
