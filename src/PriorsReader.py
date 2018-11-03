from RandomParameter import RandomParameter
from RandomParameterList import RandomParameterList
from utils import clean_tag
from lxml import etree

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
