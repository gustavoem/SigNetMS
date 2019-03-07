import sys
sys.path.insert (0, '..')

import io
import unittest
import numpy as np
from marginal_likelihood.utils import safe_exp
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
