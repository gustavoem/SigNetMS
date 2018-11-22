import sys
sys.path.insert (0, '../src/')


import os
import unittest
from ExperimentReader import read_txt_experiment_file
from ExperimentReader import read_data_experiment_file
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

    def test_read_data_experiment (self):
        """ Tests if the module can read a data experiment file. """
        data = read_data_experiment_file ("input/goodwin3.data")
        self.assertEqual (len (data), 2)
        exp0_data = data[0]
        self.assertEqual (len (exp0_data.times), 80)
        self.assertEqual (exp0_data.measure_expression, "x1")
        self.assertEqual (len (exp0_data.values), 80)
        
    def test_read_multiple_measurements (self):
        """ Tests if the module can read a data experiment with multiple
            measurements. """
        data = read_data_experiment_file ("input/goodwin3.data")
        self.assertEqual (len (data), 2)

