import sys
sys.path.insert(0, '../src/')

import unittest
import math
import numpy as np
from ODES import ODES


class TestODESMethods (unittest.TestCase):
    
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
        self.assertListEqual (list (y[0]), [100, 250.5])


    def test_integration (self):
        """ Tests if the system can be solved numerically. """
        odes = ODES ()
        # dx1 (t)/dt = x1 (t)
        # Solution is x1 (t) = exp (t)
        odes.add_equation ("x1", "x1")
        odes.define_initial_value ("x1", 1.0)
        t = np.linspace (0, 2, 11)
        y = odes.evaluate_on (t)
        for i in range (len (t)):
            analytic = math.exp (t[i])
            assert (abs (y[i] - analytic) < 1e-2)

        odes = ODES ()
        # dx1 (t)/dt = x1 (t)/2 + x2(t)/2
        # dx2 (t)/dt = x1 (t)/2 + x2(t)/2
        # Solution is x1 (t) = exp (t), x2 (t) = exp (t)
        odes.add_equation ("x1", "x1/2 + x2/2")
        odes.add_equation ("x2", "x1/2 + x2/2")
        odes.define_initial_value ("x1", 1.0)
        odes.define_initial_value ("x2", 1.0)
        y = odes.evaluate_on (t)
        for i in range (len (t)):
            analytic = math.exp (t[i])
            assert (abs (y[i][0] - analytic) < 1)
            assert (abs (y[i][1] - analytic) < 1)


    def test_equation_parameters (self):
        """ Tests if the system can store equation parameters correctly.
        """
        odes = ODES ()
        odes.add_equation ("x1", "k * x1")
        odes.define_initial_value ("x1", 1.0)
        odes.define_parameter ("k", 2)
        # Solution is x1 (t) = exp (2t)
        t = np.linspace (0, 2, 11)
        y = odes.evaluate_on (t)
        for i in range (len (t)):
            analytic = math.exp (t[i])
            assert (abs (y[i] - analytic) < 1e-2)


if __name__ == '__main__':
    unittest.main ()
