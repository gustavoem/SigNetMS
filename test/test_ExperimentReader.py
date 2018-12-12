import sys
sys.path.insert (0, '../src/')

import unittest
from ExperimentReader import read_txt_experiment_file
from Experiment import Experiment

class TestExperimentReader (unittest.TestCase):

    def test_read_txt_experiment (self):
        """ Tests if the module can read a txt experiment """
        data = read_txt_experiment_file ("input/ERK_data.txt", "ERK")
        self.assertEqual (len (data), 25)
        exp0_data = data[0]
        self.assertEqual (len (exp0_data.times), 6)
        self.assertEqual (exp0_data.measure_expression, "ERK")
        self.assertEqual (len (exp0_data.values), 6)
