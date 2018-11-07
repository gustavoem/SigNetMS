from RandomParameter import RandomParameter
from RandomParameterList import RandomParameterList
from utils import clean_tag
from lxml import etree
import sys

def read_priors_file (filename):
    tree = etree.parse (filename)
    root = tree.getroot ()

    if clean_tag (root) != "priors":
        print ("Wrong experiment data syntax. Root tag should be " + 
            "<dataset>")
        return None
    
    priors = RandomParameterList ()
    for children in root.getchildren ():
        attribs = children.attrib
        if clean_tag (children) == "prior":
            name = attribs["name"]
            a = float (attribs["a"])
            b = float (attribs["b"])
            priors.append (RandomParameter (name, a, b))
        elif clean_tag (children) == "experimental_error":
            name = attribs["name"]
            a = float (attribs["a"])
            b = float (attribs["b"])
            priors.set_experimental_error (RandomParameter (name, a ,b))
        else:
            print ("Warning: unindentified prior definition on " + filename)

    return priors


def define_sbml_params_priors (sbml, filename):
    """ Reads the priors of all sbml parameters defined in filename. """
    default_priors = read_priors_file (filename)
    priors = RandomParameterList ()

    for param in sbml.get_all_param ():
        original_name = sbml.get_original_param_name (param)
        param_prior = None
        for prior_p in default_priors.get_model_parameters ():
            if original_name == prior_p.name:
                param_prior = RandomParameter (param, prior_p.get_a (), 
                        prior_p.get_b ())
                break
        
        if param_prior == None:
            sys.exit ("ERROR! Could not find a prior for " + \
                    original_name + " in priors defined in " + filename)
        
        priors.append (param_prior)
    sigma_p = default_priors.get_experimental_error_distribution ()
    priors.set_experimental_error (sigma_p) 
    return priors
