import sys
sys.path.insert (0, '..')

import unittest
import numpy as np
from parallel_map import parallel_map

class TestParallelMap (unittest.TestCase):

    def costy_function (self, x, p1, p2, p3, p4, p5):
        y = x * p1 + p2 * p4
        for _ in range (100000):
            y = np.sqrt (np.exp (-y)) + abs (p5) * p3
        return y


    def test_parallel_implementation (self):
        """ Tests if we can compute the costy function in parallel. """
        input_nums = [np.random.randint (0, 100) for _ in range (20)]
        fun = lambda x : self.costy_function (x, 1, 20, .003, 3, -5)
        result = parallel_map (fun, input_nums, 4)
        self.assertEqual (len (result), len (input_nums))
        assert all ([r is not None for r in result])

