from Experiment import Experiment

def read_experiment_file (file_name, var):
    """ This function can read experiment files. The experiment data is 
        organized in a list of Experiment objects. """
    f = open (file_name)
    experiment_values = []
    times = []
    experiments = []
    for line in f:
        l = line.split ()
        # remove comments
        if '#' in l:
            l = l[0:l.index ('#')]
        if len (l) == 0:
            # create experiment objects
            for exp_val in experiment_values:
                exp = Experiment (times, exp_val, var)
                experiments.append (exp)
                experiment_values = []
                times = []
            continue
        
        times.append (l[0])
        if (experiment_values == []):
            nof_exp_per_line = len (l) - 1
            experiment_values = [[] for i in range (len (l) - 1)]
        for i in range (1, len (l)):
            experiment_values[i - 1].append (l[i])

    for exp_val in experiment_values:
        exp = Experiment (times, exp_val, var)
        experiments.append (exp)
        experiment_values = []
        times = []

    f.close ()
    return experiments
