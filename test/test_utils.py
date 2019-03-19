import sys
sys.path.insert (0, '..')

import io
import unittest
import numpy as np
from utils import safe_log
from utils import safe_exp
from utils import safe_exp_ratio
from utils import safe_pow_exp_ratio
from subprocess import Popen
from contextlib import redirect_stdout


class TestUtils (unittest.TestCase):
    
    def test_safe_log (self):
        """ Tests safe_log. """
        e = np.exp (1)
        self.assertEqual (safe_log (e), np.log (e))
        self.assertEqual (safe_log (1e-400), float ("-inf"))
        self.assertEqual (safe_log (0), float ("-inf"))


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


    def test_safe_pow_exp_ratio (self):
        """ Tests pow_exp_ratio function. """
        a = -1000
        b = -2000
        t = 1e-6
        r = safe_pow_exp_ratio (a, b, t)
        assert (abs (r - 0.99) < 1e-1)

