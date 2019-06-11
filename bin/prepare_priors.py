import sys
sys.path.insert (0, '..')

from model.SBML import SBML
import argparse
from lxml import etree

parser = argparse.ArgumentParser ()
parser.add_argument ("model", help="SBML file with model definition")
parser.add_argument ("output_file", help="A file to output the priors "\
        + "sketch file.")
args = parser.parse_args ()

sbml_file = args.model
output_file = args.output_file

sbml = SBML ()
sbml.load_file (sbml_file)

original_names = [sbml.get_original_param_name (p) for p in \
        sbml.get_all_param ()]
priors_names = list (set (original_names))

priors_elm = etree.Element ("priors")
for param_name in priors_names:
    param_elm = etree.SubElement (priors_elm, "prior")
    param_elm.set ("distribution", "?")
    param_elm.set ("a", "?")
    param_elm.set ("b", "?")
    param_elm.set ("name", param_name)
exp_error_elm = etree.SubElement (priors_elm, "experimental_error")
exp_error_elm.set ("distribution", "?")
exp_error_elm.set ("a", "?")
exp_error_elm.set ("b", "?")
exp_error_elm.set ("name", "Noise")
tree = etree.ElementTree (priors_elm)
tree.write (output_file, pretty_print=True, encoding='utf-8', \
        standalone=True, xml_declaration=True)
