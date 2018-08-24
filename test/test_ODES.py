import sys
sys.path.insert(0, '../src/')

import unittest
import math
import numpy as np
from ODES import ODES


class TestSBMLMethods (unittest.TestCase):
    
    def test_initial_status (self):
        """ Tests if the system starts at the initial state. """
        odes = ODES ()
        # dx1(t)/dt = x1 (t)
        # dx2(t)/dt = x2 (t)
        odes.add_equation ("x1", "x1")
        odes.add_equation ("x2", "x2")
        odes.define_initial_value ("x1", 100.0)
        odes.define_initial_value ("x2", 250.5)
        y = odes.evaluate_on ([0])
        self.assertEqual (y, [100, 250.5])


    def test_integration (self):
        """ Tests if the system can be solved numerically. """
        odes = ODES ()
        # dx1 (t)/dt = x1 (t)
        odes.add_equation ("x1", "x1")
        odes.define_initial_value ("x1", 1.0)
        # Solution is x1 (t) = exp {t^2 / 2}
        t = np.linspace (0, 2, 11)
        y = odes.evaluate_on (t)
        for i in range (len (t)):
            analytic = math.exp (t[i] ** 2 / 2)
            assert (abs (y[i] - analytic) < 1e-2)


if __name__ == '__main__':
    unittest.main ()
