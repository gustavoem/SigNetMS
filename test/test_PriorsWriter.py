import sys
sys.path.insert (0, '..')

import unittest
import filecmp
from model.PriorsReader import read_priors_file
from model.PriorsWriter import write_priors_file

class TestPriorsWriter (unittest.TestCase):

    def test_write_priors (self):
        """ Tests if it is possible to write a priors file. """
        priors = read_priors_file ('input/simple_enzymatic.priors')
        write_priors_file ('input/simple_enzymatic_out.priors', priors)
        assert (filecmp.cmp ('input/simple_enzymatic_out.priors', \
                'input/simple_enzymatic.priors'))
