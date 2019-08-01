# Prints all chemical reactions defined in an SBML file
import sys
import os
current_path = os.path.abspath (__file__)
sys.path.insert (0, '/'.join (current_path.split ('/')[:-2]))
import argparse
from model.SBML import SBML
from model.SBMLtoODES import sbml_to_odes


parser = argparse.ArgumentParser ()
parser.add_argument ("model", help="SBML file with model definition.")
args = parser.parse_args ()

sbml_file = args.model
sbml = SBML ()
sbml.load_file (sbml_file)
odes = sbml_to_odes (sbml)
odes.print_equations ()
