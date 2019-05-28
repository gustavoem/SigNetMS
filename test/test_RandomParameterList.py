import sys
sys.path.insert (0, '..')

import unittest
from model.RandomParameterList import RandomParameterList
from model.RandomParameter import RandomParameter
from distributions.Gamma import Gamma


class TestRandomParameterList (unittest.TestCase):

    def test_append (self):
        """ Tests if one can append a parameter to the list. """
        p = RandomParameter ('', Gamma (2, 2))
        theta = RandomParameterList ()
        theta.append (p)
        theta.append (p)
        theta.append (p)
        self.assertEqual (theta.get_size (), 3)


    def test_get_copy (self):
        """ Tests if an object can produce a copy of itself. """
        p1 = RandomParameter ('p1', Gamma (2, 2))
        p2 = RandomParameter ('p2', Gamma (2, 2))
        theta = RandomParameterList ()
        theta.append (p1)
        theta.append (p2)
        copy = theta.get_copy ()
        self.assertEqual (copy.get_size (), 2)
        for p in copy:
            assert (p.name == "p1" or p.name == "p2")

        p1.name = "new_name"
        for p in copy:
            assert (p.name != "new_name")


    def test_iterator (self):
        """ Tests if we can iterate through parameters. """
        p1 = RandomParameter ('p1', Gamma (2, 2))
        p2 = RandomParameter ('p2', Gamma (2, 2))
        theta = RandomParameterList ()
        theta.append (p1)
        theta.append (p2)
        for p in theta:
            assert (p.name == "p1" or p.name == "p2")


    def test_get_values (self):
        """ Tests if we can get only the values of the parameters. """
        p1 = RandomParameter ('p1', Gamma (2, 2))
        p2 = RandomParameter ('p2', Gamma (3, 2))
        p1.value = 1
        p2.value = 2
        theta = RandomParameterList ()
        theta.append (p1)
        theta.append (p2)
        values = theta.get_values ()
        self.assertEqual (values[0], 1)
        self.assertEqual (values[1], 2)
    

    def test_idexing (self):
        """ Tests if one can access a parameter with an index. """
        p1 = RandomParameter ('p1', Gamma (2, 2))
        p2 = RandomParameter ('p2', Gamma (3, 2))
        p1.value = 1
        p2.value = 2
        theta = RandomParameterList ()
        theta.append (p1)
        theta.append (p2)
        self.assertEqual (theta[0].name, 'p1')
        self.assertEqual (theta[1].name, 'p2')
        
    def test_get_experimental_error (self):
        """ Tests if we can have an experimental error on the list of
            random parameters. """
        p1 = RandomParameter ('p1', Gamma (2, 2))
        p2 = RandomParameter ('p2', Gamma (3, 2))
        p1.value = 1
        p2.value = 2
        sigma = RandomParameter ('sigma', Gamma (2, .1))
        sigma.value = 123

        theta = RandomParameterList ()
        theta.append (p1)
        theta.set_experimental_error (sigma)
        theta.append (p2)
        
        self.assertEqual (theta.get_experimental_error (), 123)
