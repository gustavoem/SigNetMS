from marginal_likelihood.RandomParameter import RandomParameter
from marginal_likelihood.RandomParameterList import RandomParameterList
from distributions.Gamma import Gamma
from distributions.Lognormal import Lognormal
from utils import clean_tag
from lxml import etree
import sys

def __create_distribution (dist_type, args):
    """ Creates a distribution dist_type with arguments args. 
        You can choose any of those options for dist_type and args,
        respectively: 
        - gamma, [k, theta]
        - lognormal, [mu, sigma] """
    if dist_type.lower () == 'gamma':
        dist = Gamma (args[0], args[1])
    elif dist_type.lower () == 'lognormal':
        dist = Lognormal (args[0], args[1])
    else:
        raise ValueError ("The specified prior distribution, " + 
                dist_type + ", is not available.")
    return dist


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
            dist_type = attribs["distribution"]
            dist = __create_distribution (dist_type, [a, b])
            priors.append (RandomParameter (name, dist))
        elif clean_tag (children) == "experimental_error":
            name = attribs["name"]
            a = float (attribs["a"])
            b = float (attribs["b"])
            gamma = Gamma (a, b)
            priors.set_experimental_error (RandomParameter (name, 
                gamma))
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
                param_prior = prior_p.copy ()
                break
        
        if param_prior == None:
            sys.exit ("ERROR! Could not find a prior for " + \
                    original_name + " in priors defined in " + filename)
        
        priors.append (param_prior)
    sigma_p = default_priors.get_experimental_error_distribution ()
    priors.set_experimental_error (sigma_p) 
    return priors
