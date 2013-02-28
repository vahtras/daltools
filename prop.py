#!/usr/bin/env python
"""Module for reading data fro property integral file AOPROPER"""

def read(n, prop="OVERLAP", propfile="AOPROPER"):
    """Read property integral"""
    import numpy
    from daltools import one
    from util import unformatted, full
    ufile = unformatted.FortranBinary(propfile)
    ufile.find(prop)
    line = ufile.data
    stars, data, symtype, label = (
          line[:8], line[8:16], line[16:24], line[24:32]
          )
    ufile.readrec()
    shape = (n, n)
    if symtype == "SQUARE  ":
        square = 1
        matsize = n*n
        buffer_ = ufile.readbuf(matsize, 'd')
        mat = full.matrix(shape)
    else:
        square = 0
        matsize = n*(n+1)/2
        buffer_ = ufile.readbuf(matsize, 'd')
        mat = numpy.array(buffer_).view(full.triangular)
        mat.anti = (symtype == "ANTISYMM")
    return mat

if __name__ == "__main__":
    import os, sys, sirifc, getopt
    from util.timing import timing
    t_all = timing("main")
    unpack = False
    verbose = False
    tmpdir = '.'
    try:
        opt, arg = getopt.getopt(
            sys.argv[1:],'t:uv',['tmpdir', 'unpack','verbose']
            )
        for o, v in opt:
            if o in ('-u','--unpack'):
                unpack = True
            if o in ('-v','--verbose'):
                verbose = True
            if o in ('-t','--tmpdir'):
                tmpdir = v
        label = arg[0]
    except IndexError:
        print "Usage: %s property" % sys.argv[0]
        sys.exit(1)
    t_read = timing("prop.read")
    propfile = os.path.join(tmpdir, 'AOPROPER')
    ifcfile = os.path.join(tmpdir, 'SIRIFC')
    n = sirifc.sirifc(name=ifcfile).nbast
    propint = read(n, prop=label, propfile=propfile)
    print t_read
    if unpack:
        t_unpack = timing("unpack")
        propsq = propint.unpack()
        print t_unpack
    if verbose:
        if unpack:
            print propsq
        else:
            print propint
    print t_all
