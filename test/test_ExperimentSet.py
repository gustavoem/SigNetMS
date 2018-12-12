import sys
sys.path.insert (0, '../src/')

import os
import unittest
from ExperimentSet import ExperimentSet
from Experiment import Experiment

class TestExperimentSet (unittest.TestCase):
    
    def test_save_to_file (self):
        """ Tests if one can save a set of experiments to a file. """ 
        times = [1, 2, 3, 4, 5]
        measure1 = 'ax + b'
        measure2 = 'e = mc ** 2'
        values1 = [123, 44, 3.2, 4, 30]
        values2 = [321, 44, 2.3, 4, 3]
        experiment1 = Experiment (times, values1, measure1)
        experiment2 = Experiment (times, values2, measure2)

        experiment_set = ExperimentSet ()
        experiment_set.add (experiment1)
        experiment_set.add (experiment2)
    
        out_file = 'tmp_exp_set.xml'
        experiment_set.save_to_file (out_file)
        
        f = open (out_file, 'r')
        data = f.read ()
        measures = [measure1, measure2]
        values = values1 + values2
        assert all ([str (t) in data for t in times])
        assert all ([m in data for m in measures])
        assert all ([str (v) in data for v in values])
