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
        self.assertListEqual (y["x1"], [100])
        self.assertListEqual (y["x2"], [250.5])


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
            assert (abs (y["x1"][i] - analytic) < 1e-2)

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
            assert (abs (y["x1"][i] - analytic) < 1)
            assert (abs (y["x2"][i] - analytic) < 1)


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
            analytic = math.exp (2 * t[i])
            assert (abs (y["x1"][i] - analytic) < 1e-2)


    def test_get_parameters (self):
        """ Tests if the system can return all parameters. """
        odes = ODES ()
        odes.add_equation ("x1", "k * x1")
        odes.define_parameter ("k", 2)
        odes.add_equation ("x2", "k2 * x2")
        odes.define_parameter ("k2", 4)
        params = odes.get_all_parameters ()
        self.assertEqual (params["k"], 2)
        self.assertEqual (params["k2"], 4)

    
    def test_get_system_jacobian (self):
        """ Tests if we can produce the jacobian of the system. """
        odes = ODES ()
        # dx1 (t)/dt = x1 (t)/2 + x2(t)/2
        # dx2 (t)/dt = x1 (t)/2 + x2(t)/2
        # Solution is x1 (t) = exp (t), x2 (t) = exp (t)
        odes.add_equation ("x1", "x1/2 + x2/2")
        odes.add_equation ("x2", "x1/2 + x2/2")
        jac_f = odes.get_system_jacobian ()
        jac = jac_f ([0, 0], [0])
        self.assertEqual (jac[0][0], .5)
        self.assertEqual (jac[0][1], .5)
        self.assertEqual (jac[1][0], .5)
        self.assertEqual (jac[1][1], .5)



if __name__ == '__main__':
    unittest.main ()
