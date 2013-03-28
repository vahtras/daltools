#!/usr/bin/env python
"""Module for reading data fro property integral file AOPROPER"""
import math
import numpy as np

def read(*args, **kwargs):
    """Read property integral"""
    from daltools import one
    from util import unformatted, full
    propfile = kwargs.get("filename", "AOPROPER")
    unpack = kwargs.get("unpack", False)
    AOPROPER = unformatted.FortranBinary(propfile)
    mat = {}
    for rec in AOPROPER:
        buf = rec.read(32, 'c')
        #stars, data, symtype, label = (
        #      buf[:8], buf[8:16], buf[16:24], buf[24:32]
        #      )
        stars = "".join(buf[:8])
        symtype = "".join(buf[16:24])
        label = "".join(buf[24:32]).strip()
        if stars == "********":
            print label, label in args
        if label.strip() in args:
            rec = next(AOPROPER)
            buffer_ = rec.read(rec.reclen/8, 'd')
            if symtype == "SQUARE  ":
                n = int(round(math.sqrt(matsize)))
                assert n*n == matsize
                mat[label] = full.init(buffer_).reshape((n, n))
            else:
                mat[label] = np.array(buffer_).view(full.triangular)
                mat[label].anti = (symtype == "ANTISYMM")
    if unpack:
        return tuple([mat[lab].unpack() for lab in args])
    else:
        return tuple([mat[lab] for lab in args])

if __name__ == "__main__":
    import os, sys, getopt
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
        proplabel = arg[0]
    except IndexError:
        print "Usage: %s property" % sys.argv[0]
        sys.exit(1)
    t_read = timing("prop.read")
    aoproper = os.path.join(tmpdir, 'AOPROPER')
    propint = read(prop=proplabel, propfile=aoproper)
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
