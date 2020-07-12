from lxml import etree

def write_priors_file (filename, priors):
    """ Saves a RandomParameterList as a priors file.
    
    Parameters
        filename: a string with the name of the priors file.
    """
    priors_elm = etree.Element ("priors")
    for param in priors.get_model_parameters ():
        param_elm = etree.SubElement (priors_elm, "prior")
        distribution = param.get_distribution ()
        param_dist_name = distribution.__class__.__name__
        param_elm.set ("distribution", param_dist_name)
        param_elm.set ("a", str (distribution.get_a ()))
        param_elm.set ("b", str (distribution.get_b ()))
        param_elm.set ("name", param.name)
    
    exp_error = priors.get_experimental_error_distribution ()
    exp_error_dist = exp_error.get_distribution ()
    exp_error_elm = etree.SubElement (priors_elm, "experimental_error")
    exp_error_dist_name = exp_error_dist.__class__.__name__
    exp_error_elm.set ("distribution", exp_error_dist_name)
    exp_error_elm.set ("a", str (exp_error_dist.get_a ()))
    exp_error_elm.set ("b", str (exp_error_dist.get_b ()))
    exp_error_elm.set ("name", "Noise")
    tree = etree.ElementTree (priors_elm)
    tree.write (filename, pretty_print=True, encoding="utf-8", \
            standalone=True, xml_declaration=True)
