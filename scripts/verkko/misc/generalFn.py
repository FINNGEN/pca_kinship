#! /usr/bin/env python
import math
import numpy as np
import re
import htmlentitydefs

##################################################################
#                   Some GlobalFunction                        #
##################################################################
'''
Length of an iterator
'''
ilen = lambda it: sum(1 for dummy in it)


##################################################################
#                   Some General Function                        #
##################################################################


def fileLength(fname):
    """
    Count the number of lines in a file
    Input : filename
    """
    i = -1
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def unique(a):
    """
    Return list with duplicate elements removed
    Input : list
    """
    return list(set(a))


def intersect(a, b):
    """
    Return the intersection of two lists
    Input : (lists,lists)
    """
    return list(set(a) & set(b))


def union(a, b):
    """
    Return the union of two lists
    Input : (lists,lists)
    """
    return list(set(a) | set(b))


def difference(a, b):
    """
    Show whats in list a which isn't in list b
    Input : (lists,lists)
    """
    return list(set(a).difference(set(b)))


def gcd(a, b):
    """
    Greatest common divisor of two numbers
    Input : (int, int)
    """
    while b:
        a, b = b, a % b
    return a


def lcm(a, b):
    """
    Lowest common factor of two numbers
    Input : (int, int)
    """
    return a*b/gcd(a, b)


def pairs(lst):
    """
    Iterate over pairs in a list (circular fashion)
    Works on any non-empty sequence, no indexing required
    Input: (list)
    Return:
     Iterator
    Usage:
    >>> list(pairs(range(5)))
    >>> [(0, 1), (1, 2), (2, 3), (3, 4), (4,0)]
    """
    i = iter(lst)
    first = prev = item = i.next()
    for item in i:
        yield prev, item
        prev = item
    yield item, first

##################################################################
#                   Some Special Function                        #
##################################################################


def cumProb(x):
    """Cumulative distribution based on rank
       Input: x (dataset)
       Output: x_cum - x variable
               y_cum - y variable
    """
    x = sorted(x, reverse=True)  # Sort the variable
    len_x = len(x)               # Length of x
    y = range(len_x-1, 0, -1)    # Range of data

    y_cum = [1]                            # First entry
    x_cum = [x[len_x-1]]                   # First entry
    # Removing duplicate entries
    for i in y:
        if x[i-1] != x[i]:
            y_cum.append(i/float(len_x))    # Non-duplicate entries
            x_cum.append(x[i-1])            # Non-duplicate entries

    # Return the values for cumulative plots
    return[x_cum, y_cum]


def sortedDictValues(adict):
    '''
    Dictionary sort by value - the simplest approach:
    Output:
    List of sorted value
    '''
    keys = adict.keys()
    keys.sort()
    return [dict[key] for key in keys]


def fxrange(limit1, limit2=None, increment=1.):
    """
    The returned value is an iterator.  Use list(fxrange) for a list.
    Range function that accepts floats (and integers).

    Usage:
    frange(-2, 2, 0.1)
    frange(10)
    frange(10, increment = 0.5)

    The returned value is an iterator.  Use list(frange) for a list.
    """

    if limit2 is None:
        limit2, limit1 = limit1, 0.
    else:
        limit1 = float(limit1)

    count = int(math.ceil((limit2 - limit1)/increment))
    return (limit1 + n*increment for n in xrange(count))


def frange(limit1, limit2=None, increment=1.):
    """
    Range function that accepts floats (and integers).

    Usage:
    frange(-2, 2, 0.1)
    frange(10)
    frange(10, increment = 0.5)
    The returned value is a list.
    """

    if limit2 is None:
        limit2, limit1 = limit1, 0.
    else:
        limit1 = float(limit1)

    count = int(math.ceil((limit2 - limit1)/increment))
    i0 = (limit1 + n*increment for n in xrange(count))
    return ["%g" % x for x in i0]


def unescape(text):
    """
    Removes HTML or XML character references and entities from a text string.
    @param text The HTML (or XML) source text.
    @return The plain text, as a Unicode string, if necessary.
    """
    def fixup(m):
        text = m.group(0)
        if text[:2] == "&#":
            # character reference
            try:
                if text[:3] == "&#x":
                    return unichr(int(text[3:-1], 16))
                else:
                    return unichr(int(text[2:-1]))
            except ValueError:
                pass
        else:
            # named entity
            try:
                text = unichr(htmlentitydefs.name2codepoint[text[1:-1]])
            except KeyError:
                pass
        return text  # leave as is
    return re.sub("&#?\w+;", fixup, text)
