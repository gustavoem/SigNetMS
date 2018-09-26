import sys
sys.path.insert (0, '../src/')

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
        self.assertEqual (exp0_data.var, "ERK")
        self.assertEqual (len (exp0_data.values), 6)

    def test_read_data_experiment (self):
        """ Tests if the module can read a data experiment file. """
        data = read_data_experiment_file ("input/goodwin3.data", "x1")
        self.assertEqual (len (data), 1)
        exp0_data = data[0]
        self.assertEqual (len (exp0_data.times), 80)
        self.asserEqual (exp0_data.var, "x1")
        self.assertEqual (len (exp0_data.values), 80)
        
