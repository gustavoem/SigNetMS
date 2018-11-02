from lxml import etree

def safe_power (a, b):
    """ Calculates a ** b safely."""
    try:
        x = a ** b
    except OverflowError:
        x = float ("inf")
    return x


def clean_tag (xmlnode):
    """ Removes the namespace of an element tag. """
    return etree.QName (xmlnode.tag).localname

