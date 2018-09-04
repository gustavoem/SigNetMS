import sys
sys.path.insert (0, '../src/')

import unittest
from ExperimentReader import read_experiment_file
from Experiment import Experiment

class TestExperimentReader (unittest.TestCase):

    def test_read_experiment (self):
        """ Tests if the module can read an experiment """
        data = read_experiment_file ('input/ERK_data.txt', 'ERK')
        self.assertEqual (len (data), 25)
        exp0_data = data[0]
        self.assertEqual (len (exp0_data.times), 6)
        self.assertEqual (exp0_data.var, "ERK")
        self.assertEqual (len (exp0_data.values), 6)
