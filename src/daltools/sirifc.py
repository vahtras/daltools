#!/usr/bin/env python
"""Module for reading SIRIUS interface file (SIRIFC)"""
import numpy
from fortran_binary import FortranBinary
from util import blocked, full


class LabelNotFound(Exception):
    pass


def get_intcode(reclen, n_floats, n_ints):
    float_size = 8
    int_size = (reclen - 4 * float_size) // n_ints
    int_codes = {4: "i", 8: "q"}
    return int_codes[int_size]


class SirIfc:
    """Read data from dalton interface file"""

    ifclabel = "SIR IPH "

    def __init__(self, name="SIRIFC"):
        self.name = name
        self.file = FortranBinary(name)
        self._cmo = None
        self._dv = None
        self._pv = None
        self._fock = None
        self._fc = None
        self._fv = None
        self._orbdiag = None

        if not self.file.find(self.ifclabel):
            raise RuntimeError(
                f"Label {self.ifclabel} not found on {name}"
            )

        rec = next(self.file)

        # Integer size from first record, 4 floats, 5 ints
        #    1) POTNUC,EMY,EACTIV,EMCSCF,ISTATE,ISPIN,NACTEL,LSYM,MS2

        try:
            self.INT = get_intcode(len(rec), 4, 5)
        except KeyError:  # pragma no cover
            # Allow Dalton 2013 version
            self.INT = get_intcode(len(rec), 4, 4)

        self.FLOAT = "d"

        self.potnuc, self.emy, self.eactive, self.emcscf = self.file.readbuf(
            4, self.FLOAT
        )
        self.istate, self.ispin, self.nactel, self.lsym = \
            self.file.readbuf(4, self.INT)

        next(self.file)
        self.nisht, self.nasht, self.nocct, self.norbt, self.nbast,\
            self.nconf,  self.nwopt, self.nwoph, self.ncdets,\
            self.ncmot, self.nnashx, self.nnashy, self.nnorbt,\
            self.n2orbt, self.nsym = \
            self.file.readbuf(15, self.INT)

        self.muld2h = numpy.array(
            self.file.readbuf(64, self.INT)
            ).reshape((8, 8))

        self.nrhf = numpy.array(self.file.readbuf(8, self.INT))
        self.nfro = numpy.array(self.file.readbuf(8, self.INT))
        self.nish = numpy.array(self.file.readbuf(8, self.INT))
        self.nash = numpy.array(self.file.readbuf(8, self.INT))
        self.norb = numpy.array(self.file.readbuf(8, self.INT))
        self.nbas = numpy.array(self.file.readbuf(8, self.INT))

        self.nelmn1, self.nelmx1, self.nelmn3, self.nelmx3, self.mctype = \
            self.file.readbuf(5, self.INT)

        self.nas1 = numpy.array(self.file.readbuf(8, self.INT))
        self.nas2 = numpy.array(self.file.readbuf(8, self.INT))
        self.nas3 = numpy.array(self.file.readbuf(8, self.INT))

        self.file.close()

    @property
    def cmo(self):
        """Read MO coefficients"""
        if self._cmo is None:
            with FortranBinary(self.name) as fb:
                fb.find(self.ifclabel)
                for _ in range(3):
                    next(fb)
                ncmot4 = max(self.ncmot, 4)
                dbl = fb.readbuf(ncmot4, self.FLOAT)
            n = 0
            self._cmo = blocked.BlockDiagonalMatrix(self.nbas, self.norb)
            for isym in range(8):
                for mo in range(self.norb[isym]):
                    for ao in range(self.nbas[isym]):
                        self._cmo.subblock[isym][ao, mo] = dbl[n]
                        n += 1
            assert n == self.ncmot
        return self._cmo

    @property
    def dv(self):
        """Get active density matrix"""
        if self._dv is None:
            with FortranBinary(self.name) as fb:
                fb.find(self.ifclabel)
                for _ in range(5):
                    next(fb)
                # mmashx = max(self.nnashx, 4)
                dbl = fb.readbuf(self.nnashx, self.FLOAT)
            self._dv = full.triangular.init(dbl)
        return self._dv

    @property
    def pv(self):
        """Get two-electron density"""
        if self._pv is None:
            with FortranBinary(self.name) as fb:
                fb.find(self.ifclabel)
                for i in range(7):
                    next(fb)
                m2ashy = max(self.nnashx ** 2, 4)
                dbl = fb.readbuf(m2ashy, self.FLOAT)
            self._pv = full.matrix((self.nnashx, self.nnashx))
            n = 0
            for i in range(self.nnashx):
                for j in range(self.nnashx):
                    self._pv[j, i] = dbl[n]
                    n += 1
            assert n == self.nnashx ** 2
        return self._pv

    @property
    def fock(self):
        """Read Fock matrix (MO)"""
        if self._fock is None:
            with FortranBinary(self.name) as fb:
                fb.find(self.ifclabel)
                for i in range(6):
                    next(fb)
                m2orbt = max(self.n2orbt, 4)
                dbl = fb.readbuf(m2orbt, self.FLOAT)
            self._fock = blocked.BlockDiagonalMatrix(self.norb, self.norb)
            n = 0
            for isym in range(8):
                for i in range(self.norb[isym]):
                    for j in range(self.norb[isym]):
                        self._fock.subblock[isym][j, i] = dbl[n]
                        n += 1
            assert n == self.n2orbt
        return self._fock

    @property
    def fc(self):
        """Read inactive Fock matrix (MO)"""
        if self._fc is None:
            with FortranBinary(self.name) as fb:
                fb.find(self.ifclabel)
                for i in range(8):
                    next(fb)
                mmorbt = max(self.nnorbt, 4)
                dbl = fb.readbuf(mmorbt, self.FLOAT)

            self._fc = blocked.triangular(self.norb)
            n = 0
            for isym in range(8):
                ij = 0
                for i in range(self.norb[isym]):
                    for j in range(i + 1):
                        self._fc.subblock[isym][i, j] = dbl[ij]
                        ij += 1
                n += ij
            assert n == self.nnorbt
        return self._fc

    @property
    def fv(self):
        """Get active Fock matrix"""
        if self._fv is None:
            with FortranBinary(self.name) as fb:
                fb.find(self.ifclabel)
                for i in range(9):
                    next(fb)
                mmorbt = max(self.nnorbt, 4)
                dbl = fb.readbuf(mmorbt, self.FLOAT)
            self._fv = blocked.triangular(self.norb)
            n = 0
            for isym in range(8):
                ij = 0
                for i in range(self.norb[isym]):
                    for j in range(i + 1):
                        self._fv.subblock[isym][i, j] = dbl[ij]
                        ij += 1
                n += ij
            assert n == self.nnorbt
        return self._fv

    @property
    def orbdiag(self):
        """Get orbital Hessian diagonal"""
        if self._orbdiag is None:
            with FortranBinary(self.name) as fb:
                if fb.find("ORBDIAG"):
                    rec = next(fb)
                    self._orbdiag = rec.read(self.nwopt, self.FLOAT)
                else:
                    raise LabelNotFound("ORBDIAG not found")
        return numpy.array(self._orbdiag)

    def __str__(self):
        retstr = ""
        retstr += f"Nuclear Potential Energy: {self.potnuc:12.6f}\n"
        retstr += f"Electronic energy       : {self.emy:12.6f}\n"
        retstr += f"Active energy           : {self.eactive:12.6f}\n"
        retstr += f"MCSCF energy            : {self.emcscf:12.6f}\n"
        retstr += f"State                   : {self.istate:d}\n"
        retstr += f"Spin                    : {self.ispin:d}\n"
        retstr += f"Active electrons        : {self.nactel:d}\n"
        retstr += f"Symmetry                : {self.lsym:d}\n"
        retstr += f"NISHT                   : {self.nisht:d}\n"
        retstr += f"NASHT                   : {self.nasht:d}\n"
        retstr += f"NOCCT                   : {self.nocct:d}\n"
        retstr += f"NORBT                   : {self.norbt:d}\n"
        retstr += f"NBAST                   : {self.nbast:d}\n"
        retstr += f"NCONF                   : {self.nconf:d}\n"
        retstr += f"NWOPT                   : {self.nwopt:d}\n"
        retstr += f"NWOPH                   : {self.nwoph:d}\n"
        retstr += f"NCDETS                  : {self.ncdets:d}\n"
        retstr += f"NCMOT                   : {self.ncmot:d}\n"
        retstr += f"NNASHX                  : {self.nnashx:d}\n"
        retstr += f"NNASHY                  : {self.nnashy:d}\n"
        retstr += f"NNORBT                  : {self.nnorbt:d}\n"
        retstr += f"N2ORBT                  : {self.n2orbt:d}\n"
        retstr += f"NSYM                    : {self.nsym:d}\n"
        retstr += "MULD2H:\n"
        for i in range(8):
            retstr += "   "
            for j in range(8):
                retstr += f" {self.muld2h[i, j]:d}"
            retstr += "\n"

        def strvec(lab, vec):
            """ String representation of vector"""
            locstr = lab + ":"
            for i in range(8):
                locstr += f" {vec[i]:d}"
            return locstr + "\n"

        retstr += strvec("NRHF", self.nrhf)
        retstr += strvec("NFRO", self.nfro)
        retstr += strvec("NISH", self.nish)
        retstr += strvec("NASH", self.nash)
        retstr += strvec("NORB", self.norb)
        retstr += strvec("NBAS", self.nbas)
        retstr += f"NELMN1                  : {self.nelmn1:d}\n"
        retstr += f"NELMX1                  : {self.nelmx1:d}\n"
        retstr += f"NELMN3                  : {self.nelmn3:d}\n"
        retstr += f"NELMX3                  : {self.nelmx3:d}\n"
        retstr += f"MCTYPE                  : {self.mctype:d}\n"
        retstr += strvec("NAS1", self.nas1)
        retstr += strvec("NAS2", self.nas2)
        retstr += strvec("NAS3", self.nas3)
        retstr += "CMO" + str(self.cmo)
        retstr += "DV\n" + str(self.dv) + "\n"
        retstr += "FOCK" + str(self.fock)
        retstr += "PV\n" + str(self.pv) + "\n"
        retstr += "FC" + str(self.fc)
        retstr += "FV" + str(self.fv)
        return retstr

    def xindx(self):
        from itertools import combinations

        ms2 = self.ispin - 1
        na = (self.nactel + ms2) // 2
        nb = (self.nactel - ms2) // 2
        astrings = tuple(combinations(range(self.nasht)[::-1], na))[::-1]
        bstrings = tuple(combinations(range(self.nasht)[::-1], nb))[::-1]
        return (
            (adet[::-1], bdet[::-1]) for bdet in bstrings for adet in astrings
        )


# For back compatibility
sirifc = SirIfc

if __name__ == "__main__":  # pragma no cover
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("sirifc", help="Sirius interface file (SIRIFC)")
    args = parser.parse_args()
    print(SirIfc(args.sirifc))
