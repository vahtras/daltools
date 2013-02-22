#!/usr/bin/env python
"""
Utility routines to generate density matrices
"""
from math import sqrt
from daltools import one, sirifc
from util import full, blocked

def h1diag(nisht, nasht):
    """Generate density from diagonalized one-electron Hamiltonian"""
    h1 = one.read('ONEHAMIL', 'AOONEINT').unpack().unblock()
    return fdiag(h1, nisht, nasht)

def fdiag(F, nisht, nasht):
    """Generate density from diagonalizing input Fock matrix"""
    C = cmo(F)
    return C2D(C, nisht, nasht)

def cmo(F, S=None, filename='AOONEINT'):
    """Return MO coefficients from diagonalization of a provided Fock matrix F,
    an optional overlap matrix S. The returned MOs are normalized and ordered 
    with respect to eigenvalue"""
    nbast = F.shape[0]
    if S is not None:
        e, V = (F/S).eigvec()
    else:
        e, V = (F).eigvec()
        S = one.read('OVERLAP', filename).unpack().unblock()
    C = full.matrix((nbast, nbast))
    N2 = (V.T*S*V).diagonal()
    for i in range(nbast):
        Ni = 1.0 / sqrt(N2[i])
        C[:, i] = V[:, i]*Ni
    #
    # Finish off with Gram Schmidt because degenerate orbitals are not
    # properly orthonormalized
    #
    C = C.GS(S)
    return C

def C2D(C, nisht, nasht):
    """Given orbitals and occupancy inactive/active 
    return inactive/active densities"""
    Ci = C[:, :nisht]
    Ca = C[:, nisht:nisht+nasht]
    Di = 2*Ci*Ci.T
    Da = Ca*Ca.T
    return Di, Da

def C1D(C, n):
    """Given orbitals C and number of occupied n return density C*C.T"""
    Ci = C[:, :n]
    Di = Ci*Ci.T
    return Di

def ifc(filename='SIRIFC', ifc_=None):
    """
    Return inactive/active densities from orbitals on interface file SIRIFC
    """
    if ifc_ is None: 
        ifc_ = sirifc.sirifc(filename)
    Di = blocked.matrix(ifc_.nbas, ifc_.nbas)
    for  isym in range(ifc_.nsym):
        if ifc_.nish[isym]:
            Ci = ifc_.cmo.subblock[isym]
            Cocci = Ci[:, :ifc_.nish[isym]]
            Di.subblock[isym] = 2*Cocci*Cocci.T
    Di = Di.unblock()
    if ifc_.nasht:
        cmoa = blocked.matrix(ifc_.nbas, ifc_.nash)
        for isym in range(ifc_.nsym):
            if ifc_.nash[isym]:
                nishi = ifc_.nish[isym]
                nocci = ifc_.nish[isym] + ifc_.nash[isym]
                cmoa.subblock[isym] = ifc_.cmo.subblock[isym][:, nishi:nocci]
        cmoa = cmoa.unblock()
        Dv = cmoa*ifc_.dv.unpack()*cmoa.T
    else:
        Dv = full.matrix((ifc_.nbast, ifc_.nbast))
    return Di, Dv

if __name__ == "__main__":
    pass
