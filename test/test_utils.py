import sys
sys.path.insert (0, '..')

import io
import unittest
import numpy as np
from marginal_likelihood.utils import safe_exp
from marginal_likelihood.utils import safe_exp_ratio
from subprocess import Popen
from contextlib import redirect_stdout


class TestUtils (unittest.TestCase):

    def test_safe_exp (self):
        """ Tests safe_exp. """
        a = safe_exp (0)
        self.assertEqual (int (a), 1)
        
        max_float = sys.float_info.max
        log_max = np.log (max_float)
        a = safe_exp (log_max)
        self.assertLess (a, float ("inf"))
        
        tmp_io = io.StringIO ()
        with redirect_stdout (tmp_io):
            a = safe_exp (log_max * 10)
        assert ('RuntimeWarning' not in tmp_io.getvalue ())


    def test_safe_exp_ratio (self):
        """ Tests the safe_exp_ratio function. """
        expected = np.exp (4 - 1.32)
        self.assertEqual (safe_exp_ratio (4, 1.32), expected)

        # if the second part is -inf then we shoulg get +inf
        expected = float ("+inf")
        self.assertEqual (safe_exp_ratio (2, float ("-inf")), expected)
