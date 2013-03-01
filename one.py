#!/usr/bin/env python
"""Module for reading data from Dalton one-electron integral file AOONEINT"""
import sys 
import struct
import numpy as np
from util import unformatted, full, blocked

def readhead(filename="AOONEINT"):
    """Read data in header of AOONEINT"""
    aooneint = unformatted.FortranBinary(filename)
    aooneint.next()
    if len(aooneint.data) == 144:
        #
        # Newer versions: title in own record
        #
        title = aooneint.data
        #
        # Next record contains MAXREP...
        #
        aooneint.next()
    else:
        #
        # Older versions: title ans MAXREP in same record
        #
        title = "".join(aooneint.readbuf(192,'c'))
    #
    # Buffer should now contains MAXREP data
    #
    nsym = aooneint.readbuf(1,'i')[0]
    #print "nsym",nsym
    naos = aooneint.readbuf(nsym,'i')
    #print "naos",naos
    potnuc = aooneint.readbuf(1,'d')[0]
    #print "naos",naos
    unlabeled = {
        "ttitle":title,
        "nsym":nsym,
        "naos":naos,
        "potnuc":potnuc
        }
    aooneint.close()
    return unlabeled

def readisordk(filename="AOONEINT"):
    """Read data under label ISORDK in AOONEINT"""
    aooneint = unformatted.FortranBinary(filename)
    table1 = aooneint.find("ISORDK")
    aooneint.next() #dummy
    aooneint.next()
    sizeofi = struct.calcsize('i')
    sizeofd = struct.calcsize('d')
    mxcent_ = (len(aooneint.data)-sizeofi)/(4*sizeofd)
    chrn_ = aooneint.readbuf(mxcent_,'d')
    nucdep = aooneint.readbuf(1,'i')[0]
    cooo_ = aooneint.readbuf(3*mxcent_,'d')
    isordk_ = {
        "table":table1,
        "chrn":chrn_,
        "nucdep":nucdep,
        "cooo":cooo_,
        }
    aooneint.close()
    return isordk_

def readscfinp(filename="AOONEINT"):
    """Read data under labl SCFINP in AOONEINT"""
    aooneint = unformatted.FortranBinary(filename)
    table2 = aooneint.find("SCFINP")
    #print table2
    scfinp_ = {}
    scfinp_["table"] = table2
    aooneint.next()
    if len(aooneint.data) == 144:
        title = aooneint.data
        aooneint.next()
    else:
        title = aooneint.readbuf(192,'c')
    scfinp_["ttitle"] = title
    nsym = aooneint.readbuf(1,'i')[0]
    scfinp_["nsym"] = nsym
    scfinp_["naos"] = aooneint.readbuf(nsym,'i')
    scfinp_["potnuc"] = aooneint.readbuf(1,'d')[0]
    kmax = aooneint.readbuf(1,'i')[0]
    scfinp_["kmax"] = kmax
    scfinp_["ncent"] = aooneint.readbuf(kmax,'i')
    nbasis = aooneint.readbuf(1,'i')[0]
    scfinp_["nbasis"] = nbasis
    scfinp_["jtran"] = aooneint.readbuf(nbasis,'i')
    scfinp_["itran"] = aooneint.readbuf(8*nbasis,'i')
    scfinp_["ctran"] = aooneint.readbuf(8*nbasis,'d')
    scfinp_["nbasis"] = aooneint.readbuf(1,'i')[0]
    scfinp_["inamn"] = aooneint.readbuf(nbasis,'i')
    scfinp_["iptyp"] = aooneint.readbuf(nbasis,'i')
    scfinp_["dpnuc"] = aooneint.readbuf(3,'d')
    nucdep = aooneint.readbuf(1,'i')[0]
    scfinp_["nucdep"] = nucdep
    scfinp_["cooo"] = aooneint.readbuf(3*nucdep,'d')
    scfinp_["ifxyz"] = aooneint.readbuf(3,'i')
    scfinp_["dummy"] = aooneint.readbuf(1,'d')[0]
    scfinp_["qpol"] = aooneint.readbuf(6,'d')
    scfinp_["qq"] = aooneint.readbuf(3,'d')
    scfinp_["jfxyz"] = aooneint.readbuf(3,'i')
    #print "scfinp_",scfinp; sys.exit(0)
    aooneint.close()
    return scfinp_

def read(label="OVERLAP", filename="AOONEINT"):
    """Read integral for label"""
    lbuf = 600
    #
    # Initialize
    #
    unlabeled = readhead(filename)
    nsym = unlabeled["nsym"]
    nbas = unlabeled["naos"]
    nnbast = 0
    for nbasi in nbas:
        nnbast += nbasi*(nbasi+1)/2
    s = full.matrix((nnbast,))
    #
    # Open file, locate label
    #
    aooneint = unformatted.FortranBinary(filename)
    labinfo = aooneint.find(label)
    #
    # Loop over records
    #
    while True:
        try:
            aooneint.next()
            buf = aooneint.readbuf(lbuf,'d')
            ibuf = aooneint.readbuf(lbuf,'i')
            length = aooneint.readbuf(1,'i')[0]
            if length < 0: raise Exception("EOD")
            for i, b in zip(ibuf[:length], buf[:length]):
                s[i-1] = b
        except Exception, inst:
            if inst[0] == "EOD":
                break
            else:
                raise
    _S = blocked.triangular(nbas)
    off = 0
    for isym in range(nsym):
        nbasi = nbas[isym]*(nbas[isym]+1)/2
        _S.subblock[isym] = np.array(s[off:off+nbasi]).view(full.triangular)
    #aooneint.close()
    return _S


if __name__ == "__main__":
    import getopt
    from util.timing import timing
    gethead = False
    getisordk = False
    getscfinp = False
    verbose = False
    unpack = False
    try:
        opt, arg = getopt.getopt(
            sys.argv[1:],'hisvu',['head','isordk','scfinp','verbose','unpack']
            )
        for o, v in opt:
            if o in ('-h','--head'):
                gethead = True
            if o in ('-i','--isordk'):
                getisordk = True
            if o in ('-s','--scfinp'):
                getscfinp = True
            if o in ('-v','--verbose'):
                verbose = True
            if o in ('-u','--unpack'):
                unpack = True
    except IndexError:
        print "Usage: %s integral" % sys.argv[0]
        sys.exit(1)

    if gethead: 
        t = timing('head')
        head = readhead()
        print "Header on AOONEINT"
        for k in head:
            print k, head[k]
        print t

    if getisordk:
        t = timing('getisrordk')
        isordk = readisordk()
        print "isordk table", isordk["table"]
        n = isordk["nucdep"]
        chrn = isordk["chrn"]
        cooo = isordk["cooo"]
        mxcent = len(chrn)
        print "nucdep=%i" % n, "mxcent=%i" % mxcent
        print full.init(chrn)[:n]
        print full.init(cooo).reshape((3, mxcent))[:, :n]
        print t

    if getscfinp:
        t = timing('scfinp')
        scfinp = readscfinp()
        for k in scfinp:
            print k, scfinp[k]
        print t

    for lab in arg:
        tread = timing("read")
        s1 = read(lab)
        print tread
        if unpack:
            t = timing("unpack")
            s2 = s1.unpack()
            print t
            t = timing("unblock")
            S = s2.unblock()
            print t
        else:
            S = s1
        if verbose:
            print lab, S
