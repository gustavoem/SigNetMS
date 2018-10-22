import sys
from sys import argv
sys.path.insert (0, '../../src/')

from lxml import etree
from Experiment import Experiment


if (len (argv) != 4):
    print ("Usage: " + argv[0] + " INPUT_FILE VAR " +
        " OUTPUT_FILE_PREFIX. Where, \n\n " + 
        "\tINPUT_FILE: the path to the input file;\n " +
        "\tVAR: the model variable measured on the experiments;\n " +
        "\tOUTPUT_FILE_PREFIX: the prefix used on the output file " +
        "names.")
else:
    filename = argv[1]
    var = argv[2]
    output_prefix = argv[3]

    f = open (filename)
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
            experiment_values[i - 1].append (float (l[i]))


    if experiment_values != []:
        for exp_val in experiment_values:
            exp = Experiment (times, exp_val, var)
            experiments.append (exp)

    f.close ()
    i = 0
    for exp in experiments:
        exp.save_to_file (output_prefix + '_' + str (i) + '.data')
        i += 1
