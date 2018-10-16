import sys
sys.path.insert (0, '../src/')

import unittest
from RandomParameterList import RandomParameterList
from RandomParameter import RandomParameter


class TestRandomParameterList (unittest.TestCase):

    def test_append (self):
        """ Tests if one can append a parameter to the list. """
        p = RandomParameter ('', 2, 2)
        theta = RandomParameterList ()
        theta.append (p)
        theta.append (p)
        theta.append (p)
        self.assertEqual (theta.get_size (), 3)


    def test_get_copy (self):
        """ Tests if an object can produce a copy of itself. """
        p1 = RandomParameter ('p1', 2, 2)
        p2 = RandomParameter ('p2', 2, 2)
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
        p1 = RandomParameter ('p1', 2, 2)
        p2 = RandomParameter ('p2', 2, 2)
        theta = RandomParameterList ()
        theta.append (p1)
        theta.append (p2)
        for p in theta:
            assert (p.name == "p1" or p.name == "p2")
