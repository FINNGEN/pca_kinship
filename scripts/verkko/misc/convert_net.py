#!/usr/bin/python
"""Convert a network file from one format to another.
Usage:
    python convertNet F_IN F_OUT {mutual} < net.in > net.out

F_IN and F_OUT are the input and output formats; see the list of
supported formats below. net.in must be in format F_IN, and the output
will be written to stdout with format F_OUT.

If the last parameter is 'mutual', only edges that go both ways will
be included.

If the same edge is encountered multiple times, the edge weight will
be the sum of all weights.

Currently supported input formats : edg
Currently supported output formats: edg, mat, gml, net (pajek)
"""

import sys
from netpython import netio

if __name__ == '__main__':
    # List of currently accepted formats.
    in_formats = ('edg',)
    out_formats = ('edg', 'mat', 'gml', 'net')

    # Get the input and output file names. If anything goes wrong,
    # print the documentation and exit.
    try:
        inputFormat = sys.argv[1]
        outputFormat = sys.argv[2]
        if len(sys.argv) > 3:
            behaviour = sys.argv[3]
        else:
            behaviour = 'simple'
        if not (inputFormat in in_formats and
                outputFormat in out_formats and 
                behaviour in ('simple', 'mutual')):
            raise
    except:
        sys.stderr.write(__doc__)
        exit(1)

    mE = (True if behaviour == 'mutual' else False)

    # Read in the network and write it out.
    netio.writeNet(netio.loadNet(sys.stdin, fileType=inputFormat, 
                                 numerical=True, mutualEdges=mE,
                                 allowSelfEdges=False),
                   sys.stdout, fileType=outputFormat)
