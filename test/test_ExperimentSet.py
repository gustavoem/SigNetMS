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
        f.close ()
        os.remove (out_file)
        measures = [measure1, measure2]
        values = values1 + values2
        assert all ([str (t) in data for t in times])
        assert all ([m in data for m in measures])
        assert all ([str (v) in data for v in values])


    def test_write_then_read (self):
        """ Tests if one can save experiments to a file and then
            read. """
        exp1 = Experiment ([1, 2, 3, 4], [.1, .2, .3, .4], 'x1')
        exp2 = Experiment ([1, 3, 9, 27], [.1, .2, .3, .4], 'x2')
        exp_set = ExperimentSet ()
        exp_set.add (exp1)
        exp_set.add (exp2)
        out_file = 'tmp_exp_set_file.xml'
        exp_set.save_to_file (out_file)
        read_exp_set = ExperimentSet (out_file)
        self.assertEqual (read_exp_set.get_size (), 2)
        exp0 = read_exp_set[0]
        exp1 = read_exp_set[1]
        self.assertListEqual (list (exp0.times), [1, 2, 3, 4])
        self.assertListEqual (list (exp0.values), [.1, .2, .3, .4])
        self.assertEqual (exp0.measure_expression, 'x1')
        self.assertListEqual (list (exp1.times), [1, 3, 9, 27])
        self.assertListEqual (list (exp1.values), [.1, .2, .3, .4])
        self.assertEqual (exp1.measure_expression, 'x2')
        os.remove (out_file)


    def test_read_experiment_set (self):
        """ Tests if the module can read a data experiment file. """
        exp_set = ExperimentSet ("input/goodwin3.data")
        self.assertEqual (exp_set.get_size (), 2)
        exp0 = exp_set[0]
        self.assertEqual (len (exp0.times), 80)
        self.assertEqual (exp0.measure_expression, "x1")
        self.assertEqual (len (exp0.values), 80)
        

    def test_read_multiple_measurements (self):
        """ Tests if the module can read a data experiment with multiple
            measurements. """
        exp_set = ExperimentSet ("input/goodwin3.data")
        self.assertEqual (exp_set[0].measure_expression, "x1")
        self.assertEqual (exp_set[1].measure_expression, "x2")
